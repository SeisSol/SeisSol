// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigBoundary.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Constants.h"
#include "Common/Marker.h"
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/configboundary.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Model/Quantities.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#include <cstdint>
#include <functional>
#include <utility>
#endif

namespace seissol::kernels {

namespace {

/// Converts `integral`, a time integral of a neighbor of the configuration `NeighborCfg`, into
/// `converted`, a time integral of the configuration `Cfg` of the cell, which the neighbor touches
/// with its side `neighborSide`; with the weights `normalStress` of the shared face if the
/// conversion reads them.
template <typename Cfg, typename NeighborCfg>
void convertNeighborIntegral(const Real<NeighborCfg>* integral,
                             Real<Cfg>* converted,
                             std::size_t neighborSide,
                             const Real<Cfg>* normalStress) {
  using CellKernels = generated::ConfigBoundaryKernels<Cfg>;
  using NeighborKernels = generated::ConfigBoundaryKernels<NeighborCfg>;

  // the constants of the kernels are in the executable on the host
  static const auto NeighborPool = seissol::Pool<NeighborCfg>::host();
  static const auto CellPool = seissol::Pool<Cfg>::host();

  alignas(Alignment) double canonical[tensor::canonicalI<NeighborCfg>::size()];

  kernel::toCanonical<NeighborCfg> toCanonical;
  toCanonical.bindGlobals(NeighborPool);
  toCanonical.I = integral;
  toCanonical.canonicalI = canonical;
  toCanonical.execute();

  // both write all of `converted`, the padding included: the neighbor kernel may read it
  if constexpr (FamilyConvertible<Cfg, NeighborCfg>) {
    // the same canonical form: both configurations are of one family
    static_assert(CellKernels::CanonicalOrder == NeighborKernels::CanonicalOrder &&
                  CellKernels::CanonicalQuantities == NeighborKernels::CanonicalQuantities);
    static_assert(tensor::canonicalI<Cfg>::size() == tensor::canonicalI<NeighborCfg>::size());

    kernel::fromCanonical<Cfg> fromCanonical;
    fromCanonical.bindGlobals(CellPool);
    fromCanonical.canonicalI = canonical;
    fromCanonical.I = converted;
    fromCanonical.execute(neighborSide);
  } else {
    // the canonical form of the family of the neighbor, which is coupled to the one of the cell
    static_assert(CellKernels::CoupledCanonicalOrder == NeighborKernels::CanonicalOrder &&
                  CellKernels::CoupledCanonicalQuantities == NeighborKernels::CanonicalQuantities);
    static_assert(tensor::coupledCanonicalI<Cfg>::size() ==
                  tensor::canonicalI<NeighborCfg>::size());

    kernel::fromCoupledCanonical<Cfg> fromCoupledCanonical;
    fromCoupledCanonical.bindGlobals(CellPool);
    fromCoupledCanonical.coupledCanonicalI = canonical;
    if constexpr (CellKernels::NormalStress) {
      assert(normalStress != nullptr);
      fromCoupledCanonical.normalStress = normalStress;
    }
    fromCoupledCanonical.I = converted;
    fromCoupledCanonical.execute(neighborSide);
  }
}

#ifdef ACL_DEVICE
/// Runs `execute` for the device kernel `krnl` of `numElements` elements with the temporary memory
/// it needs.
template <typename Kernel, typename F>
void executeWithTemporaries(Kernel& krnl,
                            std::size_t numElements,
                            seissol::parallel::runtime::StreamRuntime& runtime,
                            F&& execute) {
  krnl.numElements = numElements;
  krnl.streamPtr = runtime.stream();
  void* temporaries = nullptr;
  if constexpr (Kernel::TmpMaxMemRequiredInBytes > 0) {
    auto& device = ::device::DeviceInstance::instance();
    temporaries = device.api().allocMemAsync(Kernel::TmpMaxMemRequiredInBytes * numElements,
                                             runtime.stream());
    krnl.linearAllocator.initialize(static_cast<std::int8_t*>(temporaries));
  }
  std::invoke(std::forward<F>(execute));
  if (temporaries != nullptr) {
    ::device::DeviceInstance::instance().api().freeMemAsync(temporaries, runtime.stream());
  }
}
#endif

template <typename Cfg>
constexpr bool CanonicalOfMaterial = !generated::ConfigBoundaryKernels<Cfg>::ToCanonicalHost ||
                                     generated::ConfigBoundaryKernels<Cfg>::CanonicalQuantities ==
                                         model::MaterialOf<Cfg>::RiemannMaterial::NumQuantities;

} // namespace

template <typename Cfg>
void setNormalStressWeights(NormalStressWeights<Cfg>& weights,
                            std::size_t face,
                            const std::array<double, Cell::Dim>& normal) {
  if constexpr (generated::ConfigBoundaryKernels<Cfg>::NormalStress) {
    // the stress of the solid, in Voigt order (xx, yy, zz, xy, yz, xz)
    constexpr auto Groups = model::ElasticMaterial::PrimaryGroups;
    static_assert(model::roleKind(Groups, model::FaceRole::Traction) ==
                  model::QuantityKind::SymTensor2);
    constexpr auto Stress = model::roleOffset(Groups, model::FaceRole::Traction);
    static_assert(generated::ConfigBoundaryKernels<Cfg>::CoupledCanonicalQuantities ==
                  model::totalExtent(Groups));

    using RealT = Real<Cfg>;
    auto& faceWeights = weights.weights[face];
    std::fill(std::begin(faceWeights), std::end(faceWeights), RealT{0});
    faceWeights[Stress + 0] = static_cast<RealT>(normal[0] * normal[0]);
    faceWeights[Stress + 1] = static_cast<RealT>(normal[1] * normal[1]);
    faceWeights[Stress + 2] = static_cast<RealT>(normal[2] * normal[2]);
    faceWeights[Stress + 3] = static_cast<RealT>(2 * normal[0] * normal[1]);
    faceWeights[Stress + 4] = static_cast<RealT>(2 * normal[1] * normal[2]);
    faceWeights[Stress + 5] = static_cast<RealT>(2 * normal[0] * normal[2]);
  }
}

bool deviceConvertible(ConfigId cell, ConfigId neighbor) {
  return dispatchConfig(cell, [&](auto cellCfg) {
    return dispatchConfig(neighbor, [&](auto neighborCfg) {
      return DeviceConvertible<decltype(cellCfg), decltype(neighborCfg)>;
    });
  });
}

namespace configboundary {

recording::ConditionalKey timeKey(ConfigId neighbor, bool gts) {
  using namespace seissol::recording;
  return ConditionalKey(*KernelNames::ConfigBoundary,
                        gts ? *ComputationKind::WithGtsDerivatives
                            : *ComputationKind::WithLtsDerivatives,
                        neighbor);
}

recording::ConditionalKey toCanonicalKey(ConfigId neighbor) {
  using namespace seissol::recording;
  return ConditionalKey(*KernelNames::ConfigBoundary, *ComputationKind::None, neighbor);
}

recording::ConditionalKey fromCanonicalKey(std::size_t side) {
  using namespace seissol::recording;
  return ConditionalKey(*KernelNames::ConfigBoundary, *ComputationKind::None, *FaceId::Any, side);
}

} // namespace configboundary

template <typename Cfg>
ConfigBoundary<Cfg>::ConfigBoundary(const std::vector<ConfigId>& configs)
    : neighbors_(builtConfigCount()), devicePools_(builtConfigCount(), nullptr) {
  for (const auto config : configs) {
    if (config != configIdOf<Cfg>()) {
      neighbors_.at(config) = Neighbor{};
      empty_ = false;
    }
  }
}

template <typename Cfg>
void ConfigBoundary<Cfg>::setIntervals(double timestep,
                                       double subTimeStart,
                                       double neighborTimestep) {
  for (ConfigId config = 0; config < neighbors_.size(); ++config) {
    if (neighbors_[config].has_value()) {
      dispatchConfig(config, [&](auto neighborCfg) {
        using NeighborCfg = decltype(neighborCfg);
        const auto timeBasis = typename SolverOf<NeighborCfg>::template TimeBasis<double>(
            NeighborCfg::ConvergenceOrder);
        neighbors_[config]->timeCoeffs = timeBasis.integrate(0, timestep, timestep);
        neighbors_[config]->subtimeCoeffs =
            timeBasis.integrate(subTimeStart, timestep + subTimeStart, neighborTimestep);
      });
    }
  }
}

template <typename Cfg>
void ConfigBoundary<Cfg>::computeIntegrals(
    const CellLocalInformation& cellInformation,
    const NormalStressWeights<Cfg>& normalStress,
    const std::array<void*, Cell::NumFaces>& timeDofs,
    const std::array<real*, Cell::NumFaces>& integrationBuffer,
    std::array<real*, Cell::NumFaces>& timeIntegrated) const {
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const auto neighborConfig = cellInformation.neighborConfigIds[face];
    if (cellInformation.faceTypes[face] != FaceType::Regular ||
        neighborConfig == configIdOf<Cfg>()) {
      continue;
    }
    assert(neighborConfig < neighbors_.size() && neighbors_[neighborConfig].has_value());
    const auto& neighbor = neighbors_[neighborConfig].value();

    dispatchConfig(neighborConfig, [&](auto neighborCfg) {
      using NeighborCfg = decltype(neighborCfg);
      // a neighbor of the same configuration is skipped above
      if constexpr (!std::is_same_v<NeighborCfg, Cfg>) {
        const real* faceWeights = nullptr;
        if constexpr (generated::ConfigBoundaryKernels<Cfg>::NormalStress) {
          faceWeights = normalStress.weights[face];
        }
        computeIntegral<NeighborCfg>(
            cellInformation, face, neighbor, faceWeights, timeDofs[face], integrationBuffer[face]);
        timeIntegrated[face] = integrationBuffer[face];
      }
    });
  }
}

template <typename Cfg>
template <typename NeighborCfg>
void ConfigBoundary<Cfg>::computeIntegral(const CellLocalInformation& cellInformation,
                                          std::size_t face,
                                          const Neighbor& neighbor,
                                          const real* normalStress,
                                          const void* timeDofs,
                                          real* integrationBuffer) {
  using NeighborReal = Real<NeighborCfg>;
  if constexpr (!Convertible<Cfg, NeighborCfg>) {
    // checkConfigBoundaries admits only the faces between the configurations of one family, and
    // between a solid and a fluid
    logError() << "No conversion from the configuration"
               << configName(configValue(configIdOf<NeighborCfg>())) << "into"
               << configName(configValue(configIdOf<Cfg>())) << "was generated.";
  } else {
    alignas(Alignment) NeighborReal buffer[SolverOf<NeighborCfg>::IntegralsSize];
    const NeighborReal* integral = nullptr;
    if (cellInformation.ltsSetup.neighborBuffer(face) != BufferType::Derivatives) {
      integral = static_cast<const NeighborReal*>(timeDofs);
    } else {
      // the coefficients as in TimeCommon, in the time basis of the neighbor
      const auto& coeffs = cellInformation.ltsSetup.neighborGTSRelation(face)
                               ? neighbor.timeCoeffs
                               : neighbor.subtimeCoeffs;
      std::array<NeighborReal, NeighborCfg::ConvergenceOrder> neighborCoeffs{};
      assert(coeffs.size() == neighborCoeffs.size());
      std::transform(coeffs.begin(), coeffs.end(), neighborCoeffs.begin(), [](double coeff) {
        return static_cast<NeighborReal>(coeff);
      });
      Time<NeighborCfg> time;
      time.evaluate(neighborCoeffs.data(), static_cast<const NeighborReal*>(timeDofs), buffer);
      integral = buffer;
    }
    convertNeighborIntegral<Cfg, NeighborCfg>(
        integral, integrationBuffer, cellInformation.faceRelations[face][0], normalStress);
  }
}

template <typename Cfg>
std::size_t ConfigBoundary<Cfg>::scratchBytes(ConfigId neighbor) {
  return dispatchConfig(neighbor, [&](auto neighborCfg) -> std::size_t {
    using NeighborCfg = decltype(neighborCfg);
    if constexpr (DeviceConvertible<Cfg, NeighborCfg>) {
      using configboundary::alignScratch;
      return alignScratch(SolverOf<NeighborCfg>::IntegralsSize * sizeof(Real<NeighborCfg>)) +
             alignScratch(tensor::canonicalI<NeighborCfg>::size() * sizeof(double)) +
             alignScratch(SolverOf<Cfg>::IntegralsSize * sizeof(real));
    } else {
      return 0;
    }
  });
}

template <typename Cfg>
void ConfigBoundary<Cfg>::computeBatchedIntegrals(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) const {
#ifdef ACL_DEVICE
  using namespace seissol::recording;

  // into the canonical form, per configuration of the neighbors
  for (ConfigId config = 0; config < neighbors_.size(); ++config) {
    if (neighbors_[config].has_value()) {
      dispatchConfig(config, [&](auto neighborCfg) {
        using NeighborCfg = decltype(neighborCfg);
        if constexpr (!std::is_same_v<NeighborCfg, Cfg> && DeviceConvertible<Cfg, NeighborCfg>) {
          computeBatchedCanonical<NeighborCfg>(table, neighbors_[config].value(), runtime);
        }
      });
    }
  }

  // from the canonical form, per side of the neighbors
  if constexpr (generated::ConfigBoundaryKernels<Cfg>::Device) {
    const auto* pool = static_cast<const GlobalData<Cfg>*>(devicePools_.at(configIdOf<Cfg>()));
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      const auto key = configboundary::fromCanonicalKey(side);
      if (table.find(key) != table.end()) {
        auto& entry = table.at(key);
        auto* canonical = entry.get<double*>(inner_keys::Wp::Id::CanonicalIdofs);
        auto* integrals = entry.get<real*>(inner_keys::Wp::Id::Idofs);
        assert(pool != nullptr);
        kernel::gpu_fromCanonical<Cfg> krnl;
        krnl.bindGlobals(*pool);
        krnl.canonicalI = const_cast<const double**>(canonical->getDeviceDataPtr());
        krnl.I = integrals->getDeviceDataPtr();
        executeWithTemporaries(krnl, integrals->getSize(), runtime, [&]() { krnl.execute(side); });
      }
    }
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
template <typename NeighborCfg>
void ConfigBoundary<Cfg>::computeBatchedCanonical(
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& table,
    SEISSOL_GPU_PARAM const Neighbor& neighbor,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) const {
#ifdef ACL_DEVICE
  using namespace seissol::recording;
  using NeighborReal = Real<NeighborCfg>;

  // the time integrals of the neighbors that provide derivatives, in their time basis
  for (const bool gts : {true, false}) {
    const ConditionalKey key = configboundary::timeKey(configIdOf<NeighborCfg>(), gts);
    if (table.find(key) != table.end()) {
      auto& entry = table.at(key);
      const auto& coeffs = gts ? neighbor.timeCoeffs : neighbor.subtimeCoeffs;
      std::array<NeighborReal, NeighborCfg::ConvergenceOrder> neighborCoeffs{};
      assert(coeffs.size() == neighborCoeffs.size());
      std::transform(coeffs.begin(), coeffs.end(), neighborCoeffs.begin(), [](double coeff) {
        return static_cast<NeighborReal>(coeff);
      });
      auto* derivatives = entry.get<NeighborReal*>(inner_keys::Wp::Id::Derivatives);
      auto* integrals = entry.get<NeighborReal*>(inner_keys::Wp::Id::Idofs);
      Time<NeighborCfg> time;
      time.evaluateBatched(neighborCoeffs.data(),
                           const_cast<const NeighborReal**>(derivatives->getDeviceDataPtr()),
                           integrals->getDeviceDataPtr(),
                           integrals->getSize(),
                           runtime);
    }
  }

  const ConditionalKey key = configboundary::toCanonicalKey(configIdOf<NeighborCfg>());
  if (table.find(key) != table.end()) {
    auto& entry = table.at(key);
    auto* integrals = entry.get<NeighborReal*>(inner_keys::Wp::Id::Idofs);
    auto* canonical = entry.get<double*>(inner_keys::Wp::Id::CanonicalIdofs);
    const auto* pool =
        static_cast<const GlobalData<NeighborCfg>*>(devicePools_.at(configIdOf<NeighborCfg>()));
    assert(pool != nullptr);
    kernel::gpu_toCanonical<NeighborCfg> krnl;
    krnl.bindGlobals(*pool);
    krnl.I = const_cast<const NeighborReal**>(integrals->getDeviceDataPtr());
    krnl.canonicalI = canonical->getDeviceDataPtr();
    executeWithTemporaries(krnl, integrals->getSize(), runtime, [&]() { krnl.execute(); });
  }
#endif
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  static_assert(CanonicalOfMaterial<Cfg>);                                                         \
  template void setNormalStressWeights<Cfg>(                                                       \
      NormalStressWeights<Cfg>&, std::size_t, const std::array<double, Cell::Dim>&);               \
  template class ConfigBoundary<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels

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
#include "Common/Real.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/configboundary.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <utils/logger.h>
#include <vector>

namespace seissol::kernels {

namespace {

/// Converts `integral`, a time integral of a neighbor of the configuration `NeighborCfg`, into
/// `converted`, a time integral of the configuration `Cfg` of the cell, which the neighbor touches
/// with its side `neighborSide`.
template <typename Cfg, typename NeighborCfg>
void convertNeighborIntegral(const Real<NeighborCfg>* integral,
                             Real<Cfg>* converted,
                             std::size_t neighborSide) {
  // the same canonical form: both configurations are of one family
  static_assert(generated::ConfigBoundaryKernels<Cfg>::CanonicalOrder ==
                    generated::ConfigBoundaryKernels<NeighborCfg>::CanonicalOrder &&
                generated::ConfigBoundaryKernels<Cfg>::CanonicalQuantities ==
                    generated::ConfigBoundaryKernels<NeighborCfg>::CanonicalQuantities);
  static_assert(tensor::canonicalI<Cfg>::size() == tensor::canonicalI<NeighborCfg>::size());

  // the constants of the kernels are in the executable on the host
  static const auto NeighborPool = seissol::Pool<NeighborCfg>::host();
  static const auto CellPool = seissol::Pool<Cfg>::host();

  alignas(Alignment) double canonical[tensor::canonicalI<NeighborCfg>::size()];

  kernel::toCanonical<NeighborCfg> toCanonical;
  toCanonical.bindGlobals(NeighborPool);
  toCanonical.I = integral;
  toCanonical.canonicalI = canonical;
  toCanonical.execute();

  // writes all of `converted`, the padding included: the neighbor kernel may read it
  kernel::fromCanonical<Cfg> fromCanonical;
  fromCanonical.bindGlobals(CellPool);
  fromCanonical.canonicalI = canonical;
  fromCanonical.I = converted;
  fromCanonical.execute(neighborSide);
}

template <typename Cfg>
constexpr bool CanonicalOfMaterial = !generated::ConfigBoundaryKernels<Cfg>::Host ||
                                     generated::ConfigBoundaryKernels<Cfg>::CanonicalQuantities ==
                                         model::MaterialOf<Cfg>::RiemannMaterial::NumQuantities;

} // namespace

template <typename Cfg>
ConfigBoundary<Cfg>::ConfigBoundary(const std::vector<ConfigId>& configs)
    : neighbors_(builtConfigCount()) {
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
        computeIntegral<NeighborCfg>(
            cellInformation, face, neighbor, timeDofs[face], integrationBuffer[face]);
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
                                          const void* timeDofs,
                                          real* integrationBuffer) {
  using NeighborReal = Real<NeighborCfg>;
  if constexpr (!Convertible<Cfg, NeighborCfg>) {
    // checkConfigBoundaries admits only the faces between the configurations of one family
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
        integral, integrationBuffer, cellInformation.faceRelations[face][0]);
  }
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  static_assert(CanonicalOfMaterial<Cfg>);                                                         \
  template class ConfigBoundary<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Constants.h"
#include "Common/Real.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/configboundary.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/ConfigBoundary.h"
#include "Solver/MultipleSimulations.h"

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/Typedefs.h"
#include "Memory/GlobalData.h"
#include "Memory/MemoryAllocator.h"
#include "Parallel/Runtime/Stream.h"

#include <Device/device.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace seissol::unit_test {

namespace configboundarytest {

template <typename Cfg>
constexpr std::size_t Bases = tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg>];
template <typename Cfg>
constexpr std::size_t Quantities = tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

/// The trace of the time integral `integral` of the configuration `Cfg`, of the quantity
/// `quantity` and the simulation `sim`, on the side `Side`, in the face basis.
template <typename Cfg, std::size_t Side>
std::vector<double> trace(const Real<Cfg>* integral, std::size_t sim, std::size_t quantity) {
  const auto matrix = init::fPrT<Cfg>::template view<Side>::create(init::fPrT<Cfg>::Values[Side]);
  const auto values = init::I<Cfg>::view::create(integral);
  const bool transposed = matrix.shape(1) != Bases<Cfg>;
  const auto faceBases = transposed ? matrix.shape(1) : matrix.shape(0);
  std::vector<double> result(faceBases);
  for (std::size_t i = 0; i < matrix.shape(0); ++i) {
    for (std::size_t j = 0; j < matrix.shape(1); ++j) {
      if (matrix.isInRange(i, j)) {
        const auto face = transposed ? j : i;
        const auto basis = transposed ? i : j;
        result[face] += matrix(i, j) * multisim::multisimWrap<Cfg>(values, sim, basis, quantity);
      }
    }
  }
  return result;
}

template <typename Cfg, std::size_t... Sides>
std::vector<double> traceOnSide(const Real<Cfg>* integral,
                                std::size_t sim,
                                std::size_t quantity,
                                std::size_t side,
                                std::index_sequence<Sides...> /*sides*/) {
  std::vector<double> result;
  ((Sides == side ? (result = trace<Cfg, Sides>(integral, sim, quantity), 0) : 0), ...);
  return result;
}

/// Fills `integral`, a time integral of the configuration `Cfg`, with distinct values of order one.
template <typename Cfg>
void fillIntegral(Real<Cfg>* integral) {
  auto view = init::I<Cfg>::view::create(integral);
  for (std::size_t sim = 0; sim < Cfg::NumSimulations; ++sim) {
    for (std::size_t basis = 0; basis < Bases<Cfg>; ++basis) {
      for (std::size_t quantity = 0; quantity < Quantities<Cfg>; ++quantity) {
        // not too smooth in the basis index
        multisim::multisimWrap<Cfg>(view, sim, basis, quantity) = static_cast<Real<Cfg>>(
            std::sin(1.0 + 0.7 * static_cast<double>(basis) + 1.3 * static_cast<double>(quantity) +
                     2.1 * static_cast<double>(sim)));
      }
    }
  }
}

/// The size of the coupling matrix of a face between a cell of the configuration `Cfg` and a
/// neighbor of `NeighborCfg`, at least one.
template <typename Cfg, typename NeighborCfg>
constexpr std::size_t couplingSize() {
  if constexpr (kernels::CoupledConvertible<Cfg, NeighborCfg>) {
    return tensor::coupling<Cfg>::size(kernels::FamilyOf<NeighborCfg>);
  } else {
    return 1;
  }
}

/// The coefficients with which the quantities of a neighbor of the configuration `NeighborCfg` make
/// up the ones of the configuration `Cfg` (a row per quantity of the cell): of the same family, the
/// quantity of the same name; of the coupled family, the coupling matrix of the face, which this
/// writes to `coupling` with distinct entries of order one.
template <typename Cfg, typename NeighborCfg>
std::vector<std::vector<double>> quantityMap(Real<Cfg>* coupling) {
  std::vector<std::vector<double>> map(Quantities<Cfg>,
                                       std::vector<double>(Quantities<NeighborCfg>));
  if constexpr (kernels::CoupledConvertible<Cfg, NeighborCfg>) {
    constexpr auto Family = kernels::FamilyOf<NeighborCfg>;
    auto view = init::coupling<Cfg>::template view<Family>::create(coupling);
    // the quantities of the Riemann problem of the neighbor, in the quantities of the cell
    REQUIRE(view.shape(0) ==
            generated::ConfigBoundaryKernels<Cfg>::CoupledCanonicalQuantities[Family]);
    REQUIRE(view.shape(1) == Quantities<Cfg>);
    for (std::size_t neighborQuantity = 0; neighborQuantity < view.shape(0); ++neighborQuantity) {
      for (std::size_t quantity = 0; quantity < view.shape(1); ++quantity) {
        view(neighborQuantity, quantity) =
            static_cast<Real<Cfg>>(std::cos(0.3 + 0.9 * static_cast<double>(neighborQuantity) +
                                            1.7 * static_cast<double>(quantity)));
        map[quantity][neighborQuantity] = view(neighborQuantity, quantity);
      }
    }
  } else {
    const auto& names = model::MaterialOf<Cfg>::Quantities;
    const auto& neighborNames = model::MaterialOf<NeighborCfg>::Quantities;
    // the names of a viscoelastic material leave out its memory variables
    const auto named = std::min(Quantities<Cfg>, names.size());
    for (std::size_t quantity = 0; quantity < named; ++quantity) {
      const auto same = static_cast<std::size_t>(
          std::distance(neighborNames.begin(),
                        std::find(neighborNames.begin(), neighborNames.end(), names[quantity])));
      if (same < neighborNames.size()) {
        map[quantity][same] = 1;
      }
    }
  }
  return map;
}

/// Converts `integral`, of a neighbor of the configuration `NeighborCfg` that touches the cell with
/// its side `side`, into `converted` on the host, through `boundary`, with the coupling matrix
/// `coupling` of the shared face.
template <typename Cfg, typename NeighborCfg>
void convertOnHost(const kernels::ConfigBoundary<Cfg>& boundary,
                   Real<NeighborCfg>* integral,
                   std::size_t side,
                   Real<Cfg>* coupling,
                   Real<Cfg>* converted) {
  CellLocalInformation info{};
  info.faceTypes = {
      FaceType::Regular, FaceType::FreeSurface, FaceType::FreeSurface, FaceType::FreeSurface};
  info.neighborConfigIds = {configIdOf<NeighborCfg>(),
                            std::numeric_limits<std::uint32_t>::max(),
                            std::numeric_limits<std::uint32_t>::max(),
                            std::numeric_limits<std::uint32_t>::max()};
  info.faceRelations[0][0] = side;
  info.ltsSetup.setNeighborBuffer(0, BufferType::StepIntegrals);
  info.ltsSetup.setNeighborGTSRelation(0, true);

  const std::array<void*, Cell::NumFaces> timeDofs{integral, nullptr, nullptr, nullptr};
  const std::array<Real<Cfg>*, Cell::NumFaces> integrationBuffer{
      converted, nullptr, nullptr, nullptr};
  std::array<Real<Cfg>*, Cell::NumFaces> timeIntegrated{};
  const std::array<Real<Cfg>*, Cell::NumFaces> couplings{coupling, nullptr, nullptr, nullptr};
  boundary.computeIntegrals(info, couplings, timeDofs, integrationBuffer, timeIntegrated);
  REQUIRE(timeIntegrated[0] == converted);
}

template <typename Cfg, typename NeighborCfg>
void checkConversion() {
  using RealT = Real<Cfg>;
  using NeighborReal = Real<NeighborCfg>;

  alignas(Alignment) std::array<RealT, couplingSize<Cfg, NeighborCfg>()> coupling{};
  const auto map = quantityMap<Cfg, NeighborCfg>(coupling.data());

  // the precision of the less precise configuration
  const double tolerance =
      std::is_same_v<RealT, float> || std::is_same_v<NeighborReal, float> ? 1e-4 : 1e-10;

  kernels::ConfigBoundary<Cfg> boundary({configIdOf<Cfg>(), configIdOf<NeighborCfg>()});
  boundary.setIntervals(1.0, 0.0, 1.0);

  alignas(Alignment) std::array<NeighborReal, tensor::I<NeighborCfg>::size()> integral{};
  fillIntegral<NeighborCfg>(integral.data());

  alignas(Alignment) std::array<RealT, tensor::I<Cfg>::size()> buffer{};

  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    CAPTURE(side);
    // the conversion writes all of it
    buffer.fill(std::numeric_limits<RealT>::quiet_NaN());
    convertOnHost<Cfg, NeighborCfg>(
        boundary, integral.data(), side, coupling.data(), buffer.data());

    const auto converted = init::I<Cfg>::view::create(buffer.data());
    for (std::size_t sim = 0; sim < Cfg::NumSimulations; ++sim) {
      for (std::size_t quantity = 0; quantity < Quantities<Cfg>; ++quantity) {
        CAPTURE(quantity);
        const auto& coefficients = map[quantity];
        if (std::all_of(coefficients.begin(), coefficients.end(), [](double coefficient) {
              return coefficient == 0;
            })) {
          // nothing to take over: zero
          for (std::size_t basis = 0; basis < Bases<Cfg>; ++basis) {
            CHECK(multisim::multisimWrap<Cfg>(converted, sim, basis, quantity) == RealT{0});
          }
          continue;
        }

        // the trace on the shared face, tested with the face basis of the cell
        std::vector<double> expected;
        for (std::size_t neighborQuantity = 0; neighborQuantity < coefficients.size();
             ++neighborQuantity) {
          if (coefficients[neighborQuantity] != 0) {
            const auto part = traceOnSide<NeighborCfg>(
                integral.data(), sim, neighborQuantity, side, std::make_index_sequence<4>());
            expected.resize(std::max(expected.size(), part.size()));
            for (std::size_t face = 0; face < part.size(); ++face) {
              expected[face] += coefficients[neighborQuantity] * part[face];
            }
          }
        }
        const auto actual =
            traceOnSide<Cfg>(buffer.data(), sim, quantity, side, std::make_index_sequence<4>());
        for (std::size_t face = 0; face < actual.size(); ++face) {
          CAPTURE(face);
          const auto value = face < expected.size() ? expected[face] : 0.0;
          CHECK(actual[face] == doctest::Approx(value).epsilon(tolerance).scale(1.0));
        }
      }
    }
  }
}

#ifdef ACL_DEVICE
/// Converts a neighbor of the configuration `NeighborCfg` for all sides at once with the batched
/// kernels on the device, and compares the result with the conversion on the host, all of it.
template <typename Cfg, typename NeighborCfg>
void checkBatchedConversion() {
  using RealT = Real<Cfg>;
  using NeighborReal = Real<NeighborCfg>;
  constexpr auto Size = tensor::I<Cfg>::size();
  constexpr auto NeighborSize = tensor::I<NeighborCfg>::size();
  constexpr auto CanonicalSize = tensor::canonicalI<NeighborCfg>::size();

  // the summation orders differ between the host and the device
  const double tolerance =
      std::is_same_v<RealT, float> || std::is_same_v<NeighborReal, float> ? 1e-5 : 1e-12;

  kernels::ConfigBoundary<Cfg> boundary({configIdOf<Cfg>(), configIdOf<NeighborCfg>()});
  boundary.setIntervals(1.0, 0.0, 1.0);

  alignas(Alignment) std::array<NeighborReal, NeighborSize> integral{};
  fillIntegral<NeighborCfg>(integral.data());

  alignas(Alignment) std::array<RealT, couplingSize<Cfg, NeighborCfg>()> coupling{};
  quantityMap<Cfg, NeighborCfg>(coupling.data());

  std::vector<RealT> expected(Cell::NumFaces * Size);
  for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
    alignas(Alignment) std::array<RealT, Size> buffer{};
    convertOnHost<Cfg, NeighborCfg>(
        boundary, integral.data(), side, coupling.data(), buffer.data());
    std::copy(buffer.begin(), buffer.end(), expected.begin() + side * Size);
  }

  const auto& device = ::device::DeviceInstance::instance();
  memory::ManagedAllocator allocator;
  GlobalData<Cfg> globalData{};
  GlobalData<NeighborCfg> neighborGlobalData{};
  initializer::GlobalDataInitializerOnDevice::init<Cfg>(
      globalData, allocator, memory::Memkind::DeviceGlobalMemory);
  initializer::GlobalDataInitializerOnDevice::init<NeighborCfg>(
      neighborGlobalData, allocator, memory::Memkind::DeviceGlobalMemory);
  boundary.template setDeviceGlobalData<Cfg>(&globalData);
  boundary.template setDeviceGlobalData<NeighborCfg>(&neighborGlobalData);

  auto* deviceIntegral =
      static_cast<NeighborReal*>(device.api().allocGlobMem(NeighborSize * sizeof(NeighborReal)));
  auto* deviceCanonical =
      static_cast<double*>(device.api().allocGlobMem(CanonicalSize * sizeof(double)));
  auto* deviceConverted =
      static_cast<RealT*>(device.api().allocGlobMem(Cell::NumFaces * Size * sizeof(RealT)));
  device.api().copyTo(deviceIntegral, integral.data(), NeighborSize * sizeof(NeighborReal));
  // the coupling matrix of the shared face, as on the host
  std::vector<RealT*> deviceCoupling;
  if constexpr (kernels::CoupledConvertible<Cfg, NeighborCfg>) {
    deviceCoupling = {static_cast<RealT*>(device.api().allocGlobMem(sizeof(coupling)))};
    device.api().copyTo(deviceCoupling[0], coupling.data(), sizeof(coupling));
  }
  // the conversion writes all of it
  const std::vector<RealT> nan(Cell::NumFaces * Size, std::numeric_limits<RealT>::quiet_NaN());
  device.api().copyTo(deviceConverted, nan.data(), nan.size() * sizeof(RealT));

  {
    recording::ConditionalPointersToRealsTable table;
    std::vector<NeighborReal*> integrals{deviceIntegral};
    std::vector<double*> canonical{deviceCanonical};
    auto& toCanonical = table[kernels::configboundary::toCanonicalKey(configIdOf<NeighborCfg>())];
    toCanonical.set(recording::inner_keys::Wp::Id::Idofs, integrals);
    toCanonical.set(recording::inner_keys::Wp::Id::CanonicalIdofs, canonical);
    std::array<std::vector<RealT*>, Cell::NumFaces> converted{};
    for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
      converted[side] = {deviceConverted + side * Size};
      const auto key = kernels::FamilyConvertible<Cfg, NeighborCfg>
                           ? kernels::configboundary::fromCanonicalKey(side)
                           : kernels::configboundary::fromCoupledCanonicalKey(
                                 kernels::FamilyOf<NeighborCfg>, side);
      auto& fromCanonical = table[key];
      fromCanonical.set(recording::inner_keys::Wp::Id::CanonicalIdofs, canonical);
      fromCanonical.set(recording::inner_keys::Wp::Id::Idofs, converted[side]);
      if (!deviceCoupling.empty()) {
        fromCanonical.set(recording::inner_keys::Wp::Id::Coupling, deviceCoupling);
      }
    }

    parallel::runtime::StreamRuntime runtime;
    boundary.computeBatchedIntegrals(table, runtime);
    runtime.wait();
  }

  std::vector<RealT> actual(Cell::NumFaces * Size);
  device.api().copyFrom(actual.data(), deviceConverted, actual.size() * sizeof(RealT));
  for (std::size_t i = 0; i < actual.size(); ++i) {
    CAPTURE(i / Size);
    CAPTURE(i % Size);
    CHECK(actual[i] == doctest::Approx(expected[i]).epsilon(tolerance).scale(1.0));
  }

  for (auto* pointer : deviceCoupling) {
    device.api().freeGlobMem(pointer);
  }
  device.api().freeGlobMem(deviceConverted);
  device.api().freeGlobMem(deviceCanonical);
  device.api().freeGlobMem(deviceIntegral);
}
#endif

} // namespace configboundarytest

TEST_CASE("The canonical form of a family has the quantities of its Riemann problem") {
  forEachConfig([&](auto cfg) {
    using Cfg = decltype(cfg);
    if constexpr (generated::ConfigBoundaryKernels<Cfg>::Host) {
      CAPTURE(configName(configValue(configIdOf<Cfg>())));
      using Material = model::MaterialOf<Cfg>;
      constexpr auto Count = generated::ConfigBoundaryKernels<Cfg>::CanonicalQuantities;
      using RiemannQuantities = model::RiemannQuantitiesOf<Material>;
      REQUIRE(Count == RiemannQuantities::NumQuantities);
      REQUIRE(Count <= configboundarytest::Quantities<Cfg>);
      // the leading quantities of the configuration, in the order of the Riemann problem
      for (std::size_t quantity = 0; quantity < Count; ++quantity) {
        CHECK(Material::Quantities[quantity] == RiemannQuantities::Quantities[quantity]);
      }
    }
  });
}

TEST_CASE("A neighbor of another configuration keeps its trace on the shared face") {
  forEachConfig([&](auto cfg) {
    using Cfg = decltype(cfg);
    forEachConfig([&](auto neighborCfg) {
      using NeighborCfg = decltype(neighborCfg);
      if constexpr (!std::is_same_v<Cfg, NeighborCfg> && kernels::Convertible<Cfg, NeighborCfg>) {
        const auto cell = configName(configValue(configIdOf<Cfg>()));
        const auto neighbor = configName(configValue(configIdOf<NeighborCfg>()));
        CAPTURE(cell);
        CAPTURE(neighbor);
        configboundarytest::checkConversion<Cfg, NeighborCfg>();
      }
    });
  });
}

#ifdef ACL_DEVICE
TEST_CASE("On GPUs a neighbor of another configuration converts as on the host") {
  const auto& device = ::device::DeviceInstance::instance();
  device.api().setDevice(0);
  device.api().initialize();

  forEachConfig([&](auto cfg) {
    using Cfg = decltype(cfg);
    forEachConfig([&](auto neighborCfg) {
      using NeighborCfg = decltype(neighborCfg);
      if constexpr (!std::is_same_v<Cfg, NeighborCfg> &&
                    kernels::DeviceConvertible<Cfg, NeighborCfg>) {
        const auto cell = configName(configValue(configIdOf<Cfg>()));
        const auto neighbor = configName(configValue(configIdOf<NeighborCfg>()));
        CAPTURE(cell);
        CAPTURE(neighbor);
        configboundarytest::checkBatchedConversion<Cfg, NeighborCfg>();
      }
    });
  });

  // as in main(); otherwise, the device is torn down only by static destructors at exit
  device.api().finalize();
}
#endif

} // namespace seissol::unit_test

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
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::kernels {

namespace {

// the extents of a time integral of the configuration `Cfg`: basis functions and quantities
template <typename Cfg>
constexpr std::size_t IntegralBases = tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg>];
template <typename Cfg>
constexpr std::size_t IntegralQuantities = tensor::I<Cfg>::Shape[multisim::BasisDim<Cfg> + 1];

/// The names of the quantities of a time integral of the configuration `Cfg`; empty for the
/// quantities without a name, i.e. the memory variables of the fused layout.
template <typename Cfg>
std::vector<std::string> integralQuantities() {
  const auto& names = model::MaterialOf<Cfg>::Quantities;
  std::vector<std::string> quantities(IntegralQuantities<Cfg>);
  std::copy_n(names.begin(), std::min(names.size(), quantities.size()), quantities.begin());
  return quantities;
}

/// The trace of the volume basis of the configuration `Cfg` on the side `Side` of a cell, in the
/// face basis as the neighbor sees it: (face basis) x (volume basis).
template <typename Cfg, std::size_t Side>
Eigen::MatrixXd faceTrace() {
  const auto view = init::fPrT<Cfg>::template view<Side>::create(init::fPrT<Cfg>::Values[Side]);
  // the fused simulations store the matrix transposed
  const bool transposed = view.shape(1) != IntegralBases<Cfg>;
  const auto rows = transposed ? view.shape(1) : view.shape(0);
  Eigen::MatrixXd trace = Eigen::MatrixXd::Zero(rows, IntegralBases<Cfg>);
  for (std::size_t i = 0; i < view.shape(0); ++i) {
    for (std::size_t j = 0; j < view.shape(1); ++j) {
      if (view.isInRange(i, j)) {
        if (transposed) {
          trace(j, i) = view(i, j);
        } else {
          trace(i, j) = view(i, j);
        }
      }
    }
  }
  return trace;
}

template <typename Cfg, std::size_t... Sides>
std::array<Eigen::MatrixXd, Cell::NumFaces> faceTraces(std::index_sequence<Sides...> /*sides*/) {
  return {faceTrace<Cfg, Sides>()...};
}

template <typename Cfg>
std::array<Eigen::MatrixXd, Cell::NumFaces> faceTraces() {
  return faceTraces<Cfg>(std::make_index_sequence<Cell::NumFaces>());
}

} // namespace

NeighborConversion::NeighborConversion(ConfigId cell, ConfigId neighbor) {
  dispatchConfig(cell, [&](auto cellCfg) {
    using Cfg = decltype(cellCfg);
    dispatchConfig(neighbor, [&](auto neighborCfg) {
      using NeighborCfg = decltype(neighborCfg);

      if (Cfg::NumSimulations != NeighborCfg::NumSimulations) {
        logError() << "Neighboring cells of the configurations" << configName(configValue(cell))
                   << "and" << configName(configValue(neighbor))
                   << "fuse a different number of simulations.";
      }

      const auto names = integralQuantities<Cfg>();
      const auto neighborNames = integralQuantities<NeighborCfg>();
      quantities_.resize(names.size());
      for (std::size_t quantity = 0; quantity < names.size(); ++quantity) {
        const auto match = std::find(neighborNames.begin(), neighborNames.end(), names[quantity]);
        if (!names[quantity].empty() && match != neighborNames.end()) {
          quantities_[quantity] =
              static_cast<std::size_t>(std::distance(neighborNames.begin(), match));
        }
      }

      if (IntegralBases<NeighborCfg> > IntegralBases<Cfg>) {
        // the trace of the neighbor, projected to the face basis of the cell by truncation, is
        // matched by the minimum-norm volume polynomial of the cell
        const auto traces = faceTraces<Cfg>();
        const auto neighborTraces = faceTraces<NeighborCfg>();
        for (std::size_t side = 0; side < Cell::NumFaces; ++side) {
          const auto& trace = traces[side];
          const Eigen::MatrixXd rightInverse =
              trace.transpose() * (trace * trace.transpose()).inverse();
          const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> lift =
              rightInverse * neighborTraces[side].topRows(trace.rows());
          lift_[side].assign(lift.data(), lift.data() + lift.size());
        }
      }
    });
  });
}

template <typename Cfg, typename NeighborCfg>
void NeighborConversion::apply(const Real<NeighborCfg>* integral,
                               Real<Cfg>* converted,
                               std::size_t neighborSide) const {
  using RealT = Real<Cfg>;
  constexpr auto Bases = IntegralBases<Cfg>;
  constexpr auto NeighborBases = IntegralBases<NeighborCfg>;

  const auto from = init::I<NeighborCfg>::view::create(integral);
  // the padding included: the neighbor kernel may read it
  std::fill_n(converted, tensor::I<Cfg>::size(), RealT{0});
  auto to = init::I<Cfg>::view::create(converted);

  const auto& lift = lift_[neighborSide];
  for (std::size_t sim = 0; sim < Cfg::NumSimulations; ++sim) {
    for (std::size_t quantity = 0; quantity < quantities_.size(); ++quantity) {
      if (!quantities_[quantity].has_value()) {
        continue;
      }
      const auto neighborQuantity = quantities_[quantity].value();
      for (std::size_t basis = 0; basis < Bases; ++basis) {
        double value = 0;
        if (lift.empty()) {
          if (basis < NeighborBases) {
            value = multisim::multisimWrap<NeighborCfg>(from, sim, basis, neighborQuantity);
          }
        } else {
          for (std::size_t neighborBasis = 0; neighborBasis < NeighborBases; ++neighborBasis) {
            value +=
                lift[basis * NeighborBases + neighborBasis] *
                multisim::multisimWrap<NeighborCfg>(from, sim, neighborBasis, neighborQuantity);
          }
        }
        multisim::multisimWrap<Cfg>(to, sim, basis, quantity) = static_cast<RealT>(value);
      }
    }
  }
}

template <typename Cfg>
ConfigBoundary<Cfg>::ConfigBoundary(const std::vector<ConfigId>& configs)
    : neighbors_(builtConfigCount()) {
  for (const auto config : configs) {
    if (config != configIdOf<Cfg>()) {
      neighbors_.at(config) = Neighbor{NeighborConversion(configIdOf<Cfg>(), config), {}, {}};
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
      using NeighborReal = Real<NeighborCfg>;
      if constexpr (!std::is_same_v<NeighborCfg, Cfg>) {
        alignas(Alignment) NeighborReal buffer[SolverOf<NeighborCfg>::IntegralsSize];
        const NeighborReal* integral = nullptr;
        if (cellInformation.ltsSetup.neighborBuffer(face) != BufferType::Derivatives) {
          integral = static_cast<const NeighborReal*>(timeDofs[face]);
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
          time.evaluate(
              neighborCoeffs.data(), static_cast<const NeighborReal*>(timeDofs[face]), buffer);
          integral = buffer;
        }
        neighbor.conversion.template apply<Cfg, NeighborCfg>(
            integral, integrationBuffer[face], cellInformation.faceRelations[face][0]);
        timeIntegrated[face] = integrationBuffer[face];
      }
    });
  }
}

#define SEISSOL_INSTANTIATE(Cfg) template class ConfigBoundary<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels

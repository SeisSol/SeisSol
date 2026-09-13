// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_SOLVER_H_
#define SEISSOL_SRC_KERNELS_STP_SOLVER_H_

#include "Kernels/Common.h"

#include <cstddef>
#include <variant>

namespace seissol::numerical {
template <typename>
class LegendreBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {
class Local;
class Neighbor;
} // namespace seissol::kernels::solver::linearck

namespace seissol::tensor {
struct spaceTimePredictor;
} // namespace seissol::tensor

namespace seissol::kernels::solver::stp {

class Spacetime;
class Time;

struct STPLocalData;

struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = Time;
  using LocalKernelT = linearck::Local;
  using NeighborKernelT = linearck::Neighbor;

  template <typename RealT>
  using TimeBasis = seissol::numerical::LegendreBasis<RealT>;

  /// Whether the face flux solvers are built from the flux itself rather than
  /// from a Godunov state. A Godunov state is a state, so the solver it
  /// builds maps the state onto itself; a solver that transports more than
  /// the state cannot use one, and assembles the flux of the face normal
  /// instead.
  static constexpr bool FluxSolverFromTable = false;
  static constexpr std::size_t IntegralsSize = tensor::I::size();
  static constexpr std::size_t DerivativesSize = kernels::size<tensor::spaceTimePredictor>();

  using LocalData = STPLocalData;
  using NeighborData = std::monostate;
};

} // namespace seissol::kernels::solver::stp
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

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

  /// The basis of whatever else a solver keeps an expansion of. A solver that
  /// transports only its state expands only in the basis its recursion
  /// produces; one that transports more has to win those coefficients from
  /// samples, and a monomial basis loses five digits at order six doing it.
  template <typename RealT>
  using ExtraTimeBasis = TimeBasis<RealT>;

  /// Whether the face flux solvers are built from the flux itself rather than
  /// from a Godunov state. A Godunov state is a state, so the solver it
  /// builds maps the state onto itself; a solver that transports more than
  /// the state cannot use one, and assembles the flux of the face normal
  /// instead.
  /// Whether a cell has to keep the integrals it handed its neighbours.
  ///
  /// A face flux that is scaled with the larger of the two sides' wave speeds
  /// has to be formed where both sides are, which is the neighbouring
  /// integration -- and there a cell's own integrals are not passed in. They
  /// cannot be rebuilt from its derivatives either, because what it
  /// transports is not an expansion of its state alone. So the cell keeps
  /// them, whatever the cluster relations of its faces would otherwise ask
  /// for: that includes a cell whose faces are all boundaries, which would
  /// have no buffer at all.
  static constexpr bool RequiresOwnIntegrals = false;
  static constexpr bool FluxSolverFromTable = false;
  static constexpr std::size_t IntegralsSize = tensor::I::size();
  static constexpr std::size_t DerivativesSize = kernels::size<tensor::spaceTimePredictor>();

  using LocalData = STPLocalData;
  using NeighborData = std::monostate;
};

} // namespace seissol::kernels::solver::stp
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

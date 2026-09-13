// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_

#include "GeneratedCode/tensor.h"

#include <cstddef>
#include <variant>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {

struct LinearLocalData;

class Spacetime;
class Time;
class Local;
class Neighbor;

struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = Time;
  using LocalKernelT = Local;
  using NeighborKernelT = Neighbor;

  template <typename RealT>
  using TimeBasis = seissol::numerical::MonomialBasis<RealT>;

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
  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ>();

  using LocalData = LinearLocalData;
  using NeighborData = std::monostate;
};

} // namespace seissol::kernels::solver::linearck
#endif // SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_

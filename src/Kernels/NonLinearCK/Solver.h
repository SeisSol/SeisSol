// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_

#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"

#include <cstddef>
#include <variant>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
template <typename>
class LegendreBasis;
} // namespace seissol::numerical

namespace seissol::tensor {
struct transportDer;
}

namespace seissol::kernels::solver::nonlinearck {

struct NonLinearLocalData;
struct NonLinearNeighborData;

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

  /// The basis of whatever else a solver keeps an expansion of. A solver that
  /// transports only its state expands only in the basis its recursion
  /// produces; one that transports more has to win those coefficients from
  /// samples, and a monomial basis loses five digits at order six doing it.
  template <typename RealT>
  using ExtraTimeBasis = seissol::numerical::LegendreBasis<RealT>;

  /// A cell hands its neighbours one tensor: the time-integrated state it
  /// couples through, the time-integrated stress, and the wave speed the
  /// dissipation is scaled with. The stress is carried rather than recomputed
  /// because the flux is nonlinear in time -- the integral of the stress over
  /// a timestep is not the stress of the integrated strain once the internal
  /// variables move within the step -- and carrying it means neither side of a
  /// face ever needs the other's material. `I` is therefore wider than `Q`,
  /// and the internal variables are not part of it: they carry no flux, so no
  /// neighbour reads them.
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
  static constexpr bool RequiresOwnIntegrals = true;
  static constexpr bool FluxSolverFromTable = true;
  static constexpr std::size_t IntegralsSize = tensor::I::size();

  /// The expansion of the state, and behind it the expansion of what the cell
  /// carries beyond it. The second has no recursion to come out of, so it is
  /// projected from the samples of the step -- and a neighbour on a coarser
  /// cluster reconstructs a subinterval of the stress from it rather than
  /// rebuilding the stress out of this cell's material.
  /// Size of the expansion of what a cell carries beyond its state.
  static constexpr std::size_t DerivativesSize =
      yateto::computeFamilySize<tensor::dQ>() + kernels::familySize<tensor::transportDer>();

  using LocalData = NonLinearLocalData;
  using NeighborData = NonLinearNeighborData;
};

} // namespace seissol::kernels::solver::nonlinearck
#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_

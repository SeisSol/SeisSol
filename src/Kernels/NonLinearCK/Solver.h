// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_

#include "GeneratedCode/init.h"
#include "GeneratedCode/quantities.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/TimeCoefficients.h"
#include "Numerical/TimeBasis.h"

#include <algorithm>
#include <cstddef>
#include <variant>
#include <yateto/InitTools.h>

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

  /// The bases of every expansion this solver keeps. The state expands in the
  /// basis its recursion produces; what it carries beyond the state has to be
  /// won from samples, and a monomial basis loses five digits at order six
  /// doing it.
  template <typename RealT>
  using TimeBasis = seissol::numerical::CompoundTimeBasis<TimeCoefficients,
                                                          seissol::numerical::MonomialBasis<RealT>,
                                                          seissol::numerical::LegendreBasis<RealT>>;

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

  /// How a coarser cluster folds one step's integrals into what it already
  /// has. Every column is a time integral and therefore a sum -- except the
  /// one that carries the bound a face scales its dissipation with, which is
  /// a maximum over the step. The larger of two bounds is a bound; their sum
  /// is not one, it is merely larger, and it would grow with the cluster
  /// ratio.
  static void accumulate(real* accumulated, const real* step) {
    const auto bound = generated::TransportBoundColumn;
    const auto entry = init::I::index(0, bound);
#pragma omp simd
    for (std::size_t dof = 0; dof < IntegralsSize; ++dof) {
      accumulated[dof] += step[dof];
    }
    accumulated[entry] = std::max(accumulated[entry] - step[entry], step[entry]);
  }

  using LocalData = NonLinearLocalData;
  using NeighborData = NonLinearNeighborData;
};

} // namespace seissol::kernels::solver::nonlinearck
#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_

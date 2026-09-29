// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_SOLVER_H_

#include "GeneratedCode/quantities.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Common.h"
#include "Kernels/TimeCoefficients.h"
#include "Numerical/TimeBasis.h"

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
                                                          TimeQuadrature,
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

  /// Whether the predictor has to sample its flux over the step instead of
  /// integrating it in closed form. A solver that samples is handed the rule to
  /// sample with; forming one for a solver that does not would cost a little
  /// per step and answer nothing.
  static constexpr bool RequiresTimeQuadrature = true;

  static constexpr bool FluxSolverFromTable = true;

  /// Which boundary conditions this solver's kernels implement.
  ///
  /// A face without a neighbor needs a ghost rule that folds into the pair of
  /// flux matrices, and for a transported tensor that means a mirror of the
  /// transported quantities: outflow needs none, a free surface has one, and a
  /// fault carries the flux its friction imposes. The conditions that impose a
  /// value instead act on the quantities of the Riemann problem, which is a
  /// narrower thing than what this solver transports -- so their kernels are
  /// not generated for this layout at all (`codegen/kernels/nodalbc.py`), and
  /// `Setup.h` has no rule to build the pair from. Stating it here is what
  /// turns that into a rejection at startup with the reason, rather than an
  /// error once the setup reaches the offending face.
  static constexpr FaceTypeSupport implementsFaceType(FaceType faceType) {
    if (faceType == FaceType::FreeSurfaceGravity) {
      return faceTypeUnsupported("the surface elevation is imposed on the quantities of the "
                                 "Riemann problem, which is not what this solver transports");
    }
    if (faceType == FaceType::Dirichlet) {
      return faceTypeUnsupported("the Dirichlet datum is folded into the rows of the local flux "
                                 "solver, and a transported tensor has more than those rows");
    }
    if (faceType == FaceType::Analytical) {
      return faceTypeUnsupported("this predictor evaluates no time-dependent boundary condition");
    }
    return faceTypeSupported();
  }

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

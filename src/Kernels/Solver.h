// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_SOLVER_H_
#define SEISSOL_SRC_KERNELS_SOLVER_H_

#include "Config.h"
#include "Equations/Datastructures.h"
#include "Kernels/Data.h"
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Kernels/Precision.h"
#include "Kernels/STP/Solver.h"
#include "Kernels/TimeCoefficients.h"
#include "Numerical/TimeBasis.h"

#include <vector>

// IWYU pragma: begin_exports

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "Kernels/LinearCKAnelastic/Local.h"
#include "Kernels/LinearCKAnelastic/Neighbor.h"
#include "Kernels/LinearCKAnelastic/Time.h"
#elif defined(SEISSOL_KERNELS_NONLINEARCK)
#include "Kernels/NonLinearCK/Local.h"
#include "Kernels/NonLinearCK/Neighbor.h"
#include "Kernels/NonLinearCK/Time.h"
#elif defined(SEISSOL_KERNELS_STP)
#include "Kernels/LinearCK/Local.h"
#include "Kernels/LinearCK/Neighbor.h"
#include "Kernels/STP/Time.h"
#else
#include "Kernels/LinearCK/Local.h"
#include "Kernels/LinearCK/Neighbor.h"
#include "Kernels/LinearCK/Time.h"
#endif

// IWYU pragma: end_exports

namespace seissol::kernels {

// some typename shortcuts

using Solver = model::MaterialT::Solver;

using Time = Solver::TimeKernelT;
using Spacetime = Solver::SpacetimeKernelT;
using Local = Solver::LocalKernelT;
using Neighbor = Solver::NeighborKernelT;

inline Solver::TimeBasis<real> timeBasis() {
  return Solver::TimeBasis<real>(Config::ConvergenceOrder);
}

/// Coefficients that integrate over [start, end] of a step of `timestep`.
inline TimeCoefficients timeIntegrate(double start, double end, double timestep) {
  return timeBasis().integrate(start, end, timestep);
}

/// Coefficients that evaluate at a point of a step of `timestep`.
inline TimeCoefficients timePoint(double position, double timestep) {
  return timeBasis().point(position, timestep);
}

/// One set of coefficients per point, which is what a face interpolating in
/// time asks for. Where a flat concatenation used to be sliced by hand, the
/// points are addressed by their index.
inline std::vector<TimeCoefficients> timeCollocate(const std::vector<double>& points,
                                                   double timestep) {
  std::vector<TimeCoefficients> coefficients;
  coefficients.reserve(points.size());
  for (const auto& point : points) {
    coefficients.push_back(timePoint(point, timestep));
  }
  return coefficients;
}

/// Everything a predictor needs of one step of `timestep`, formed once for the
/// step: the coefficients that integrate over it, and -- where the solver has
/// to sample a flux it cannot integrate -- the rule it samples with, whose
/// nodes include the ends of the step so that the chain of nodes covers it
/// without a gap and the first node is the state itself.
inline TimeStepCoefficients timeStepCoefficients(double timestep) {
  TimeStepCoefficients coefficients{timeIntegrate(0, timestep, timestep), {}};
  if constexpr (Solver::RequiresTimeQuadrature) {
    coefficients.quadrature = timeBasis().quadratureWithEndpoints(timestep);
  }
  return coefficients;
}

} // namespace seissol::kernels
#endif // SEISSOL_SRC_KERNELS_SOLVER_H_

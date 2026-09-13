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

#include <algorithm>
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

inline Solver::ExtraTimeBasis<real> extraTimeBasis() {
  return Solver::ExtraTimeBasis<real>(Config::ConvergenceOrder);
}

namespace detail {
template <typename Operation>
inline TimeCoefficients bothBases(Operation&& operation) {
  TimeCoefficients coefficients{};
  const auto state = operation(timeBasis());
  const auto extra = operation(extraTimeBasis());
  std::copy(state.begin(), state.end(), coefficients.state.begin());
  std::copy(extra.begin(), extra.end(), coefficients.extra.begin());
  return coefficients;
}
} // namespace detail

/// Coefficients that integrate over [start, end] of a step of `timestep`.
inline TimeCoefficients timeIntegrate(double start, double end, double timestep) {
  return detail::bothBases(
      [&](const auto& basis) { return basis.integrate(start, end, timestep); });
}

/// Coefficients that evaluate at a point of a step of `timestep`.
inline TimeCoefficients timePoint(double position, double timestep) {
  return detail::bothBases([&](const auto& basis) { return basis.point(position, timestep); });
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

} // namespace seissol::kernels
#endif // SEISSOL_SRC_KERNELS_SOLVER_H_

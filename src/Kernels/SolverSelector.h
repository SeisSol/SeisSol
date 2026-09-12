// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_SOLVERSELECTOR_H_
#define SEISSOL_SRC_KERNELS_SOLVERSELECTOR_H_

#include "Common/Typedefs.h"
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Kernels/NonLinearCK/Solver.h"
#include "Kernels/STP/Solver.h"

namespace seissol::kernels {

/// Maps the configured solver onto its implementation. Which solvers a
/// material may be built with is checked when the build is configured; nothing
/// here restricts the choice.
template <SolverType Solver>
struct SolverSelector;

template <>
struct SolverSelector<SolverType::LinearCK> {
  using Type = solver::linearck::Solver;
};

template <>
struct SolverSelector<SolverType::LinearCKAnelastic> {
  using Type = solver::linearckanelastic::Solver;
};

template <>
struct SolverSelector<SolverType::NonLinearCK> {
  using Type = solver::nonlinearck::Solver;
};

template <>
struct SolverSelector<SolverType::STP> {
  using Type = solver::stp::Solver;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SOLVERSELECTOR_H_

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
#include "Kernels/STP/Solver.h"

namespace seissol::kernels {

/// Maps the solver of a configuration onto its implementation. Which solvers a
/// material may be built with is checked when the build is configured; nothing
/// here restricts the choice.
template <SolverType Solver, typename Cfg>
struct SolverSelector;

template <typename Cfg>
struct SolverSelector<SolverType::LinearCK, Cfg> {
  using Type = solver::linearck::Solver<Cfg>;
};

template <typename Cfg>
struct SolverSelector<SolverType::LinearCKAnelastic, Cfg> {
  using Type = solver::linearckanelastic::Solver<Cfg>;
};

template <typename Cfg>
struct SolverSelector<SolverType::STP, Cfg> {
  using Type = solver::stp::Solver<Cfg>;
};

/// The solver a configuration advances its cells with.
template <typename Cfg>
using SolverOf = typename SolverSelector<Cfg::Solver, Cfg>::Type;

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SOLVERSELECTOR_H_

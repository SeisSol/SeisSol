// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SOLVERSELECTOR_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SOLVERSELECTOR_H_

#include "Core/Common/Typedefs.h"
#include "WavePropagation/Solver/LinearCK/Solver.h"
#include "WavePropagation/Solver/LinearCKAnelastic/Solver.h"
#include "WavePropagation/Solver/STP/Solver.h"

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
using SolverOf = SolverSelector<Cfg::Solver, Cfg>::Type;

} // namespace seissol::kernels

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SOLVERSELECTOR_H_

// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SOLVER_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SOLVER_H_

#include "Core/Common/Real.h"
#include "Core/Numerical/TimeBasis.h"
#include "WavePropagation/Solver/Data.h"
#include "WavePropagation/Solver/LinearCK/Solver.h"
#include "WavePropagation/Solver/LinearCKAnelastic/Solver.h"
#include "WavePropagation/Solver/STP/Solver.h"
#include "WavePropagation/Solver/SolverSelector.h"

// IWYU pragma: begin_exports

// the kernels of the solvers the configurations of the build advance their cells with
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "WavePropagation/Solver/LinearCKAnelastic/Local.h"
#include "WavePropagation/Solver/LinearCKAnelastic/Neighbor.h"
#include "WavePropagation/Solver/LinearCKAnelastic/Time.h"
#endif
#if defined(SEISSOL_KERNELS_LINEARCK) || defined(SEISSOL_KERNELS_STP)
#include "WavePropagation/Solver/LinearCK/Local.h"
#include "WavePropagation/Solver/LinearCK/Neighbor.h"
#endif
#ifdef SEISSOL_KERNELS_LINEARCK
#include "WavePropagation/Solver/LinearCK/Time.h"
#endif
#ifdef SEISSOL_KERNELS_STP
#include "WavePropagation/Solver/STP/Time.h"
#endif

// IWYU pragma: end_exports

namespace seissol::kernels {

// the kernels of the solver of a configuration `Cfg`

template <typename Cfg>
using Time = SolverOf<Cfg>::TimeKernelT;
template <typename Cfg>
using Spacetime = SolverOf<Cfg>::SpacetimeKernelT;
template <typename Cfg>
using Local = SolverOf<Cfg>::LocalKernelT;
template <typename Cfg>
using Neighbor = SolverOf<Cfg>::NeighborKernelT;

/// The time basis the solver of the configuration `Cfg` expands in.
template <typename Cfg>
using TimeBasis = SolverOf<Cfg>::template TimeBasis<Real<Cfg>>;

template <typename Cfg>
TimeBasis<Cfg> timeBasis() {
  return TimeBasis<Cfg>(Cfg::ConvergenceOrder);
}

} // namespace seissol::kernels
#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SOLVER_H_

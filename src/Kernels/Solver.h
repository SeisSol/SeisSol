// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_SOLVER_H_
#define SEISSOL_SRC_KERNELS_SOLVER_H_

#include "Common/Real.h"
#include "Kernels/Data.h"
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Kernels/STP/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Numerical/TimeBasis.h"

// IWYU pragma: begin_exports

#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "Kernels/LinearCKAnelastic/Local.h"
#include "Kernels/LinearCKAnelastic/Neighbor.h"
#include "Kernels/LinearCKAnelastic/Time.h"
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

// the kernels of the solver of a configuration `Cfg`

template <typename Cfg>
using Time = typename SolverOf<Cfg>::TimeKernelT;
template <typename Cfg>
using Spacetime = typename SolverOf<Cfg>::SpacetimeKernelT;
template <typename Cfg>
using Local = typename SolverOf<Cfg>::LocalKernelT;
template <typename Cfg>
using Neighbor = typename SolverOf<Cfg>::NeighborKernelT;

/// The time basis the solver of the configuration `Cfg` expands in.
template <typename Cfg>
using TimeBasis = typename SolverOf<Cfg>::template TimeBasis<Real<Cfg>>;

template <typename Cfg>
TimeBasis<Cfg> timeBasis() {
  return TimeBasis<Cfg>(Cfg::ConvergenceOrder);
}

} // namespace seissol::kernels
#endif // SEISSOL_SRC_KERNELS_SOLVER_H_

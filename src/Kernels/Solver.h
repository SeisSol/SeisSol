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
#include "Kernels/LinearCK/Solver.h"
#include "Kernels/LinearCKAnelastic/Solver.h"
#include "Kernels/Precision.h"
#include "Kernels/STP/Solver.h"
#include "Numerical/TimeBasis.h"

// IWYU pragma: begin_exports

#ifdef USE_VISCOELASTIC2
#include "Kernels/LinearCKAnelastic/LocalBase.h"
#include "Kernels/LinearCKAnelastic/NeighborBase.h"
#include "Kernels/LinearCKAnelastic/TimeBase.h"
#elif defined(USE_POROELASTIC)
#include "Kernels/LinearCK/LocalBase.h"
#include "Kernels/LinearCK/NeighborBase.h"
#include "Kernels/STP/TimeBase.h"
#else
#include "Kernels/LinearCK/LocalBase.h"
#include "Kernels/LinearCK/NeighborBase.h"
#include "Kernels/LinearCK/TimeBase.h"
#endif

// IWYU pragma: end_exports

namespace seissol::kernels {

// some typename shortcuts

using Solver = typename model::MaterialT::Solver;

using Time = typename Solver::TimeKernelT;
using Spacetime = typename Solver::SpacetimeKernelT;
using Local = typename Solver::LocalKernelT;
using Neighbor = typename Solver::NeighborKernelT;

inline typename Solver::TimeBasis<real> timeBasis() {
  return Solver::TimeBasis<real>(Config::ConvergenceOrder);
}

} // namespace seissol::kernels
#endif // SEISSOL_SRC_KERNELS_SOLVER_H_

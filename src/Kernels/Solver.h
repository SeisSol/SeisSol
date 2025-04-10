#pragma once

#include <Kernels/LinearCK/Solver.h>
#include <Kernels/LinearCKAnelastic/Solver.h>
#include <Kernels/STP/Solver.h>

// IWYU pragma: begin_exports

#ifdef USE_VISCOELASTIC2
#include <Kernels/LinearCKAnelastic/LocalBase.h>
#include <Kernels/LinearCKAnelastic/NeighborBase.h>
#include <Kernels/LinearCKAnelastic/TimeBase.h>
#elif defined(USE_STP)
#include <Kernels/LinearCK/LocalBase.h>
#include <Kernels/LinearCK/NeighborBase.h>
#include <Kernels/LinearCK/TimeBase.h>
#include <Kernels/STP/TimeBase.h>
#else
#include <Kernels/LinearCK/LocalBase.h>
#include <Kernels/LinearCK/NeighborBase.h>
#include <Kernels/LinearCK/TimeBase.h>
#endif

// IWYU pragma: end_exports

namespace seissol::kernels {

#ifdef USE_VISCOELASTIC2
using Solver = solver::linearckanelastic::Solver;
#elif defined(USE_STP)
using Solver = solver::stp::Solver;
#else
using Solver = solver::linearck::Solver;
#endif

using Time = typename Solver::TimeKernelT;
using Spacetime = typename Solver::SpacetimeKernelT;
using Local = typename Solver::LocalKernelT;
using Neighbor = typename Solver::NeighborKernelT;

} // namespace seissol::kernels

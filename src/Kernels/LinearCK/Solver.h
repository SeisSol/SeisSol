#pragma once

#include <cstddef>
namespace seissol::kernels::solver::linearck {

class Spacetime;
class Time;
class Local;
class Neighbor;

class Solver {
    using SpacetimeKernelT = Spacetime;
    using TimeKernelT = Time;
    using LocalKernelT = Local;
    using NeighborKernelT = Neighbor;
};

} // namespace seissol::kernels::solver::linearck

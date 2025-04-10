#pragma once

#include <cstddef>

namespace seissol::kernels::solver::linearck {
class Time;
class Local;
class Neighbor;
} // namespace seissol::kernels::solver::linearck

namespace seissol::kernels::solver::stp {

class Spacetime;

struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = linearck::Time;
  using LocalKernelT = linearck::Local;
  using NeighborKernelT = linearck::Neighbor;
};

} // namespace seissol::kernels::solver::stp

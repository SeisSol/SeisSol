// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_SOLVER_H_
#define SEISSOL_SRC_KERNELS_STP_SOLVER_H_

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
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

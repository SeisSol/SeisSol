// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

#include <cstddef>
namespace seissol::kernels::solver::linearckanelastic {

class Spacetime;
class Time;
class Local;
class Neighbor;

struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = Time;
  using LocalKernelT = Local;
  using NeighborKernelT = Neighbor;
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

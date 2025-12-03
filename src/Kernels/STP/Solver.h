// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_SOLVER_H_
#define SEISSOL_SRC_KERNELS_STP_SOLVER_H_

#include "Kernels/Common.h"

#include <cstddef>

namespace seissol::numerical {
template <typename>
class LegendreBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {
class Local;
class Neighbor;
} // namespace seissol::kernels::solver::linearck

namespace seissol::tensor {
class spaceTimePredictor;
} // namespace seissol::tensor

namespace seissol::kernels::solver::stp {

class Spacetime;
class Time;

struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = Time;
  using LocalKernelT = linearck::Local;
  using NeighborKernelT = linearck::Neighbor;

  template <typename RealT>
  using TimeBasis = seissol::numerical::LegendreBasis<RealT>;

  static constexpr std::size_t DerivativesSize = kernels::size<tensor::spaceTimePredictor>();
};

} // namespace seissol::kernels::solver::stp
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

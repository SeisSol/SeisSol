// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_SOLVER_H_
#define SEISSOL_SRC_KERNELS_STP_SOLVER_H_

#include <Kernels/Common.h>
#include <cstddef>

namespace seissol::numerical {
template <typename>
class LegendreBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {
template <typename>
class Local;
template <typename>
class Neighbor;
} // namespace seissol::kernels::solver::linearck

namespace seissol::tensor {
template <typename>
class spaceTimePredictor;
} // namespace seissol::tensor

namespace seissol::kernels::solver::stp {

template <typename>
class Spacetime;
template <typename>
class Time;

struct Solver {
  template <typename Cfg>
  using SpacetimeKernelT = Spacetime<Cfg>;
  template <typename Cfg>
  using TimeKernelT = Time<Cfg>;
  template <typename Cfg>
  using LocalKernelT = linearck::Local<Cfg>;
  template <typename Cfg>
  using NeighborKernelT = linearck::Neighbor<Cfg>;

  template <typename RealT>
  using TimeBasis = seissol::numerical::LegendreBasis<RealT>;

  static constexpr std::size_t DerivativesSize = kernels::size<tensor::spaceTimePredictor<Cfg>>();
};

} // namespace seissol::kernels::solver::stp
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

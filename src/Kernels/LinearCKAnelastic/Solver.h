// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

#include <GeneratedCode/tensor.h>
#include <cstddef>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearckanelastic {

template <typename>
class Spacetime;
template <typename>
class Time;
template <typename>
class Local;
template <typename>
class Neighbor;

struct Solver {
  template <typename Cfg>
  using SpacetimeKernelT = Spacetime<Cfg>;

  template <typename Cfg>
  using TimeKernelT = Time<Cfg>;

  template <typename Cfg>
  using LocalKernelT = Local<Cfg>;

  template <typename Cfg>
  using NeighborKernelT = Neighbor<Cfg>;

  template <typename RealT>
  using TimeBasis = seissol::numerical::MonomialBasis<RealT>;

  template <typename Cfg>
  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ<Cfg>>();
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

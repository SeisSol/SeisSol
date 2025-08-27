// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

#include <cstddef>
#include <generated_code/tensor.h>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
} // namespace seissol::numerical

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

  template <typename RealT>
  using TimeBasis = seissol::numerical::MonomialBasis<RealT>;

  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ>();
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

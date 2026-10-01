// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"

#include <cstddef>
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

template <typename Cfg>
struct AnelasticLocalData;
template <typename Cfg>
struct AnelasticNeighborData;

/// The solver as a configuration `Cfg` runs it. Its kernels are those of the configuration of the
/// build.
template <typename Cfg>
struct Solver {
  using SpacetimeKernelT = Spacetime;
  using TimeKernelT = Time;
  using LocalKernelT = Local;
  using NeighborKernelT = Neighbor;

  template <typename RealT>
  using TimeBasis = seissol::numerical::MonomialBasis<RealT>;

  static constexpr FaceTypeSupport implementsFaceType(FaceType /*faceType*/) {
    return faceTypeSupported();
  }

  /// The memory variables keep the mechanism in a tensor dimension of their own, apart from the
  /// quantity axis.
  static constexpr bool FusedMechanisms = false;

  static constexpr std::size_t IntegralsSize = tensor::I<Cfg>::size();
  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ<Cfg>>();

  using LocalData = AnelasticLocalData<Cfg>;
  using NeighborData = AnelasticNeighborData<Cfg>;
};

} // namespace seissol::kernels::solver::linearckanelastic
#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_SOLVER_H_

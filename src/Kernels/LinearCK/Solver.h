// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"

#include <cstddef>
#include <variant>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {

template <typename Cfg>
struct LinearLocalData;

class Spacetime;
class Time;
class Local;
class Neighbor;

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

  /// The memory variables of a material with relaxation share the quantity axis with its other
  /// quantities.
  static constexpr bool FusedMechanisms = true;

  static constexpr std::size_t IntegralsSize = tensor::I<Cfg>::size();
  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ<Cfg>>();

  using LocalData = LinearLocalData<Cfg>;
  using NeighborData = std::monostate;
};

} // namespace seissol::kernels::solver::linearck
#endif // SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_

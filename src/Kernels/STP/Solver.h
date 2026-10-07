// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_SOLVER_H_
#define SEISSOL_SRC_KERNELS_STP_SOLVER_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Common.h"

#include <cstddef>
#include <variant>

namespace seissol::numerical {
template <typename>
class LegendreBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {
template <typename Cfg>
class Local;
template <typename Cfg>
class Neighbor;
} // namespace seissol::kernels::solver::linearck

namespace seissol::kernels::solver::stp {

template <typename Cfg>
class Spacetime;
template <typename Cfg>
class Time;

template <typename Cfg>
struct STPLocalData;

/// The solver as a configuration `Cfg` runs it, with the kernels of `Cfg`.
template <typename Cfg>
struct Solver {
  using SpacetimeKernelT = Spacetime<Cfg>;
  using TimeKernelT = Time<Cfg>;
  using LocalKernelT = linearck::Local<Cfg>;
  using NeighborKernelT = linearck::Neighbor<Cfg>;

  template <typename RealT>
  using TimeBasis = seissol::numerical::LegendreBasis<RealT>;

  static constexpr FaceTypeSupport implementsFaceType(FaceType faceType) {
    if (faceType == FaceType::FreeSurfaceGravity || faceType == FaceType::NonlinearDirichlet) {
      // The surface elevation is built up from the Taylor derivative family dQ, which the
      // space-time predictor does not produce; so is the inner state of the nonlinear Dirichlet
      // boundary at a time.
      return faceTypeUnsupported("the space-time predictor provides no Taylor derivatives");
    }
    return faceTypeSupported();
  }

  /// Only the linear Cauchy-Kovalevskaya solver puts memory variables on the quantity axis.
  static constexpr bool FusedMechanisms = false;

  static constexpr std::size_t IntegralsSize = tensor::I<Cfg>::size();
  static constexpr std::size_t DerivativesSize = kernels::size<tensor::spaceTimePredictor<Cfg>>();

  using LocalData = STPLocalData<Cfg>;
  using NeighborData = std::monostate;
};

} // namespace seissol::kernels::solver::stp
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

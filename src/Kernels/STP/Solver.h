// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_SOLVER_H_
#define SEISSOL_SRC_KERNELS_STP_SOLVER_H_

#include "Initializer/BasicTypedefs.h"
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
struct spaceTimePredictor;
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

  static constexpr FaceTypeSupport implementsFaceType(FaceType faceType) {
    if (faceType == FaceType::FreeSurfaceGravity) {
      // The surface elevation is built up from the Taylor derivative family dQ, which the
      // space-time predictor does not produce.
      return faceTypeUnsupported("the space-time predictor provides no Taylor derivatives");
    }
    return faceTypeSupported();
  }

  static constexpr std::size_t IntegralsSize = tensor::I::size();
  static constexpr std::size_t DerivativesSize = kernels::size<tensor::spaceTimePredictor>();
};

} // namespace seissol::kernels::solver::stp
#endif // SEISSOL_SRC_KERNELS_STP_SOLVER_H_

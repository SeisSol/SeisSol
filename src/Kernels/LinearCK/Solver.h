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
#include "Model/OperatorLayout.h"

#include <cstddef>
#include <variant>
#include <yateto/InitTools.h>

namespace seissol::numerical {
template <typename>
class MonomialBasis;
} // namespace seissol::numerical

namespace seissol::kernels::solver::linearck {

struct LinearLocalData;

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

  static constexpr FaceTypeSupport implementsFaceType(FaceType faceType) {
    if (NodalFlux &&
        (faceType == FaceType::Dirichlet || faceType == FaceType::FreeSurfaceGravity)) {
      // The map of these conditions is folded into the flux matrix of the face
      // (foldDirichlet, foldFreeSurfaceGravity), and a flux that reads the
      // material at the nodes of the face does not apply that matrix.
      return faceTypeUnsupported("the boundary map is folded into a flux matrix that a flux "
                                 "varying along the face does not read");
    }
    return faceTypeSupported();
  }

  static constexpr std::size_t IntegralsSize = tensor::I::size();
  static constexpr std::size_t DerivativesSize = yateto::computeFamilySize<tensor::dQ>();

  using LocalData = LinearLocalData;
  using NeighborData = std::monostate;
};

} // namespace seissol::kernels::solver::linearck
#endif // SEISSOL_SRC_KERNELS_LINEARCK_SOLVER_H_

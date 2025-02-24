// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GravitationalFreeSurfaceBC.h"
#include <Common/Constants.h>
#include <Initializer/BasicTypedefs.h>
#include <cstddef>
#include <generated_code/kernel.h>
#include <generated_code/tensor.h>
#include <utility>

namespace seissol {

std::pair<long long, long long>
    GravitationalFreeSurfaceBc::getFlopsDisplacementFace(unsigned int face, FaceType faceType) {
  long long hardwareFlops = 0;
  long long nonZeroFlops = 0;

  constexpr auto NumberOfNodes = nodal::tensor::nodes2D::Shape[0];

  // initialize integral of displacement
  hardwareFlops += 1 * NumberOfNodes;
  nonZeroFlops += 1 * NumberOfNodes;

  // Before adjusting the range of the loop, check range of loop in computation!
  for (std::size_t order = 1; order < ConvergenceOrder + 1; ++order) {
    constexpr auto FlopsPerQuadpoint = 4 + // Computing coefficient
                                       6 + // Updating displacement
                                       2;  // Updating integral of displacement
    constexpr auto FlopsUpdates = FlopsPerQuadpoint * NumberOfNodes;

    nonZeroFlops += kernel::projectDerivativeToNodalBoundaryRotated::nonZeroFlops(order - 1, face) +
                    FlopsUpdates;
    hardwareFlops +=
        kernel::projectDerivativeToNodalBoundaryRotated::hardwareFlops(order - 1, face) +
        FlopsUpdates;
  }

  // Two rotations: One to face-aligned, one to global
  hardwareFlops += 2 * kernel::rotateFaceDisplacement::HardwareFlops;
  nonZeroFlops += 2 * kernel::rotateFaceDisplacement::NonZeroFlops;

  return {nonZeroFlops, hardwareFlops};
}
} // namespace seissol

// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GravitationalFreeSurfaceBC.h"
#include <Common/Constants.h>
#include <GeneratedCode/kernel.h>
#include <GeneratedCode/tensor.h>
#include <Initializer/BasicTypedefs.h>
#include <Solver/MultipleSimulations.h>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace seissol {

template<typename Cfg>
std::pair<std::uint64_t, std::uint64_t>
    GravitationalFreeSurfaceBc<Cfg>::getFlopsDisplacementFace(unsigned int face, FaceType faceType) {
  std::uint64_t hardwareFlops = 0;
  std::uint64_t nonZeroFlops = 0;

  constexpr std::uint64_t NumberOfNodes =
      static_cast<std::uint64_t>(nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisFunctionDimension]) *
      multisim::NumSimulations;

  // initialize integral of displacement
  hardwareFlops += 1 * NumberOfNodes;
  nonZeroFlops += 1 * NumberOfNodes;

  // Before adjusting the range of the loop, check range of loop in computation!
  for (std::size_t order = 1; order < Cfg::ConvergenceOrder + 1; ++order) {
    constexpr auto FlopsPerQuadpoint = 4 + // Computing coefficient
                                       6 + // Updating displacement
                                       2;  // Updating integral of displacement
    constexpr auto FlopsUpdates = FlopsPerQuadpoint * NumberOfNodes;

    nonZeroFlops += kernel::projectDerivativeToNodalBoundaryRotated<Cfg>::nonZeroFlops(order - 1, face) +
                    FlopsUpdates;
    hardwareFlops +=
        kernel::projectDerivativeToNodalBoundaryRotated<Cfg>::hardwareFlops(order - 1, face) +
        FlopsUpdates;
  }

  // Two rotations: One to face-aligned, one to global
  hardwareFlops += 2 * kernel::rotateFaceDisplacement<Cfg>::HardwareFlops;
  nonZeroFlops += 2 * kernel::rotateFaceDisplacement<Cfg>::NonZeroFlops;

  return {nonZeroFlops, hardwareFlops};
}

#define _H_(cfg) template class GravitationalFreeSurfaceBc<cfg>;
#include "ConfigInclude.h"

} // namespace seissol

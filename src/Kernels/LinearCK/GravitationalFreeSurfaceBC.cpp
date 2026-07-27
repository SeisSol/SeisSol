// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GravitationalFreeSurfaceBC.h"

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Monitoring/Metric.h"
#include "Solver/MultipleSimulations.h"

#include <cstddef>
#include <cstdint>

namespace seissol {

PerformanceEstimate GravitationalFreeSurfaceBc::metrics(int8_t face, FaceType /*faceType*/) {
  PerformanceEstimate estimate;

  constexpr std::uint64_t NumberOfNodes =
      static_cast<std::uint64_t>(nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension]) *
      multisim::NumSimulations;

  // initialize integral of displacement
  estimate.hardwareFlop += 1 * NumberOfNodes;
  estimate.nonzeroFlop += 1 * NumberOfNodes;

  // Before adjusting the range of the loop, check range of loop in computation!
  for (std::size_t order = 1; order < ConvergenceOrder + 1; ++order) {
    constexpr auto FlopsPerQuadpoint = 4 + // Computing coefficient
                                       6 + // Updating displacement
                                       2;  // Updating integral of displacement
    constexpr auto FlopsUpdates = FlopsPerQuadpoint * NumberOfNodes;

    estimate += PerformanceEstimate::fromKernel<kernel::projectDerivativeToNodalBoundaryRotated>(
        order - 1, face);

    estimate.hardwareFlop += FlopsUpdates;
    estimate.nonzeroFlop += FlopsUpdates;
  }

  // Two rotations: One to face-aligned, one to global
  estimate += PerformanceEstimate::fromKernel<kernel::rotateFaceDisplacement>() * 2;

  return estimate;
}
} // namespace seissol

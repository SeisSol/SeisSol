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
#include "Solver/MultipleSimulations.h"

#include <cstddef>
#include <cstdint>
#include <utility>

namespace seissol {

std::pair<std::uint64_t, std::uint64_t>
    GravitationalFreeSurfaceBc::getFlopsDisplacementFace(unsigned int face) {
  std::uint64_t hardwareFlops = 0;
  std::uint64_t nonZeroFlops = 0;

  hardwareFlops += kernel::fsgKernel::hardwareFlops(face);
  nonZeroFlops += kernel::fsgKernel::nonZeroFlops(face);

  return {nonZeroFlops, hardwareFlops};
}
} // namespace seissol

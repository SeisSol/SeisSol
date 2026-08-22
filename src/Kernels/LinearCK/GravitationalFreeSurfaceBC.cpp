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

PerformanceEstimate GravitationalFreeSurfaceBc::metrics(int8_t face) {
  return PerformanceEstimate::fromKernel<kernel::fsgKernel>(face);
}

} // namespace seissol

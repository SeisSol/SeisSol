// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_STP_DATA_H_
#define SEISSOL_SRC_KERNELS_STP_DATA_H_

#include "Kernels/Common.h"
#include "Kernels/Precision.h"

#include <cstddef>

namespace seissol::tensor {
struct ET;
struct Zinv;
} // namespace seissol::tensor

namespace seissol::kernels::solver::stp {
// TODO: remove zeroGuard when only initialized where relevant

struct STPLocalData {
  real sourceMatrix[zeroGuard(kernels::size<tensor::ET>())]{};

  // currently hard-coded to poroelasticity
  // NOLINTNEXTLINE
  real G[13]{};

  // currently hard-coded to poroelasticity
  // NOLINTNEXTLINE
  real Zinv[zeroGuard(kernels::familySize<tensor::Zinv>())]{};

  // preferrably double; will be compared closely against the "default" timestep width almost all
  // the time
  double typicalTimeStepWidth{};
};

} // namespace seissol::kernels::solver::stp

#endif // SEISSOL_SRC_KERNELS_STP_DATA_H_

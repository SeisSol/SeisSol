// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DATA_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DATA_H_

#include "Kernels/Common.h"
#include "Kernels/Precision.h"

#include <cstddef>

namespace seissol::tensor {
struct ET;
} // namespace seissol::tensor

namespace seissol::kernels::solver::linearck {
// TODO: remove zeroLengthArrayHandler when only initialized where relevant

struct LinearLocalData {
  real sourceMatrix[zeroLengthArrayHandler(kernels::size<tensor::ET>())]{};
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DATA_H_

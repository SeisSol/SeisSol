// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DATA_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DATA_H_

#include "Common/Real.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"

#include <cstddef>

namespace seissol::kernels::solver::linearck {
// TODO: remove zeroGuard when only initialized where relevant

template <typename Cfg>
struct LinearLocalData {
  Real<Cfg> sourceMatrix[zeroGuard(kernels::size<tensor::ET<Cfg>>())]{};
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DATA_H_

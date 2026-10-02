// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_DEVICEAUX_KERNELSAUX_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_DEVICEAUX_KERNELSAUX_H_

#include "GeneratedCode/init.h"
#include "Kernels/Precision.h"

namespace seissol::kernels::time::aux {
template <typename Cfg>
void taylorSum(std::size_t count,
               Real<Cfg>** target,
               const Real<Cfg>** source,
               const Real<Cfg>* coeffs,
               void* stream);
} // namespace seissol::kernels::time::aux

#endif // SEISSOL_SRC_KERNELS_LINEARCK_DEVICEAUX_KERNELSAUX_H_

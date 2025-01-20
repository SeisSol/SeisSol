// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_TIMEBASE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_TIMEBASE_H_

#include "generated_code/kernel.h"

namespace seissol {
namespace kernels {
class TimeBase {
  protected:
  kernel::derivative m_krnlPrototype;
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_TIMEBASE_H_

// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_

#include "generated_code/kernel.h"

namespace seissol {
namespace kernels {
class NeighborBase {
  protected:
  kernel::neighbourFluxExt m_nfKrnlPrototype;
  kernel::neighbour m_nKrnlPrototype;
  dynamicRupture::kernel::nodalFlux m_drKrnlPrototype;
};
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_KERNELS_NEIGHBORBASE_H_

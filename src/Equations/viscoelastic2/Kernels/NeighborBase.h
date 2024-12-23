// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alexander Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

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

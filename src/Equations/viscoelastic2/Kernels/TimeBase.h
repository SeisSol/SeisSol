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

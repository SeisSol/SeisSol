// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_

#define DMO_INCREMENT 1
#define DMO_STREAM(IN, OUT) *(OUT) = *(IN);

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_

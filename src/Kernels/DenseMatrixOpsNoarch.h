// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_

#define DMO_INCREMENT 1
#define DMO_BROADCAST(IN, OUT) real OUT = *(IN);
#define DMO_STREAM(IN, OUT) *(OUT) = *(IN);
#define DMO_SXT(S, X, Y) *(Y) = (S) * *(X);
#define DMO_SXTYP(S, X, Y) *(Y) += (S) * *(X);
#define DMO_XYMST(S, X, Y, Z) *(Z) = (*(X)-*(Y)) * (S);
#define DMO_XYMSTZP(S, X, Y, Z) *(Z) += (*(X)-*(Y)) * (S);

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_


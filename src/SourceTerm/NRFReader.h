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

#ifndef SEISSOL_SRC_SOURCETERM_NRFREADER_H_
#define SEISSOL_SRC_SOURCETERM_NRFREADER_H_

#include "NRF.h"

namespace seissol::sourceterm {
void readNRF(const char* filename, NRF& nrf);
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_NRFREADER_H_

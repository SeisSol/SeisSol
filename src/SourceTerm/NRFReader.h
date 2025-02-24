// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOURCETERM_NRFREADER_H_
#define SEISSOL_SRC_SOURCETERM_NRFREADER_H_

#include "NRF.h"

namespace seissol::sourceterm {
void readNRF(const char* filename, NRF& nrf);
} // namespace seissol::sourceterm

#endif // SEISSOL_SRC_SOURCETERM_NRFREADER_H_

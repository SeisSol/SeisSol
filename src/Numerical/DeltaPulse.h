// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_DELTAPULSE_H_
#define SEISSOL_SRC_NUMERICAL_DELTAPULSE_H_

#include "Common/Marker.h"

namespace seissol::deltaPulse {

SEISSOL_HOSTDEVICE inline real deltaPulse(real time, real timeStep) {

  if (time > 0 && time <= timeStep) {
    return (1 / timeStep);
  } else {
    return 0;
  }
}

} // namespace seissol::deltaPulse

#endif // SEISSOL_SRC_NUMERICAL_DELTAPULSE_H_

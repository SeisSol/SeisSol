// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_NUMERICAL_DELTAPULSE_H_
#define SEISSOL_SRC_NUMERICAL_DELTAPULSE_H_

namespace seissol::deltaPulse {

inline real deltaPulse(real time, real timeStep) {

  if (time > 0 && time <= timeStep) {
    return (1 / timeStep);
  } else {
    return 0;
  }
}

} // namespace seissol::deltaPulse

#endif // SEISSOL_SRC_NUMERICAL_DELTAPULSE_H_

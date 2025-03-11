#ifndef SEISSOL_DELTAPULSE_H
#define SEISSOL_DELTAPULSE_H

namespace seissol::deltaPulse {

inline real deltaPulse(real time, real timeStep, real surfaceArea) {

  if (time > 0 && time <= timeStep*1.00001) {
    return (1 / timeStep) * (1 / surfaceArea);
  } else {
    return 0;
  }
}

} // namespace seissol::deltaPulse

#endif // SEISSOL_DELTAPULSE_H

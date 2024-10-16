#ifndef SEISSOL_DELTAPULSE_H
#define SEISSOL_DELTAPULSE_H

namespace seissol::deltaPulse {

inline real deltaPulse(real time, real timeStep, real elementSurfaceArea){

    if(time > 0 && time <= timeStep){
        return (1 / timeStep) * (1 / elementSurfaceArea);
    }
    else{
        return 0;
    }
}

} // namespace seissol::deltaPulse

#endif // SEISSOL_DELTAPULSE_H
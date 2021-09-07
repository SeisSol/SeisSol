#include "GravitationalFreeSurfaceBC.h"
#include "SeisSol.h"

namespace seissol {

double getGravitationalAcceleration() {
  return SeisSol::main.getGravitationSetup().acceleration;
}

} // namespace seissol

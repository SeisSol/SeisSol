#include "Misc.h"

#include <Geometry/MeshDefinition.h>
#include <cmath>

namespace seissol::dr::misc {
template <>
double power<8>(double base) {
  const double baseSquared = base * base;
  const double baseToTheFour = baseSquared * baseSquared;
  const double baseToTheEight = baseToTheFour * baseToTheFour;
  return baseToTheEight;
}

void computeStrikeAndDipVectors(const VrtxCoords normal, VrtxCoords strike, VrtxCoords dip) {
  // compute normalized strike vector
  const auto strikeInvLength = 1.0 / std::sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
  strike[0] = normal[1] * strikeInvLength;
  strike[1] = -normal[0] * strikeInvLength;
  strike[2] = 0.0;

  // compute normalized dip vector
  dip[0] = -1.0 * strike[1] * normal[2];
  dip[1] = strike[0] * normal[2];
  dip[2] = strike[1] * normal[0] - strike[0] * normal[1];
  const auto dipInvLength = 1.0 / std::sqrt(dip[0] * dip[0] + dip[1] * dip[1] + dip[2] * dip[2]);
  dip[0] *= dipInvLength;
  dip[1] *= dipInvLength;
  dip[2] *= dipInvLength;
}
} // namespace seissol::dr::misc
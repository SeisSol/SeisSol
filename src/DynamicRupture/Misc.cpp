#include "Misc.h"

#include <cmath>

namespace seissol::dr::misc {
template <>
real power<8>(real base) {
  real baseSquared = base * base;
  real baseToTheFour = baseSquared * baseSquared;
  real baseToTheEight = baseToTheFour * baseToTheFour;
  return baseToTheEight;
}
real magnitude(real x, real y) { return std::sqrt(x * x + y * y); }
real asinh(real x) { return std::log(x + std::sqrt(misc::power<2>(x) + 1.0)); }
} // namespace seissol::dr::misc
#include "Misc.h"

#include <cmath>

namespace seissol::dr::misc {
template <>
double power<8>(double base) {
  double baseSquared = base * base;
  double baseToTheFour = baseSquared * baseSquared;
  double baseToTheEight = baseToTheFour * baseToTheFour;
  return baseToTheEight;
}
real magnitude(real x, real y) { return std::sqrt(x * x + y * y); }
double asinh(double x) { return std::log(x + std::sqrt(x * x + 1.0)); }
} // namespace seissol::dr::misc
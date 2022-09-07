#ifndef SEISSOL_REGULARIZEDYOFFE_H
#define SEISSOL_REGULARIZEDYOFFE_H

namespace seissol::regularizedYoffe {
/**
 * Implementation of the regularized Yoffe function defined in Appendix of Tinti et al. (2005)
 */
inline real regularizedYoffe(real time, real tauS, real tauR) {
  real k = 2.0 / (M_PI * tauR * tauS * tauS);
  // c1 to c6 are analytical functions used for building the regularized Yoffe function
  auto c1 = [&]() {
    return (0.5 * time + 0.25 * tauR) * std::sqrt(time * (tauR - time)) +
           (time * tauR - tauR * tauR) * std::asin(std::sqrt(time / tauR)) -
           0.75 * tauR * tauR * std::atan(std::sqrt((tauR - time) / time));
  };

  auto c2 = [&] { return 0.375 * M_PI * tauR * tauR; };

  auto c3 = [&]() {
    return (tauS - time - 0.5 * tauR) * std::sqrt((time - tauS) * (tauR - time + tauS)) +
           tauR * (2 * tauR - 2 * time + 2 * tauS) * std::asin(std::sqrt((time - tauS) / tauR)) +
           1.5 * tauR * tauR * std::atan(std::sqrt((tauR - time + tauS) / (time - tauS)));
  };

  auto c4 = [&]() {
    // 2 typos fixed in the second term compared with Tinti et al. 2005
    return (-tauS + 0.5 * time + 0.25 * tauR) *
               std::sqrt((time - 2.0 * tauS) * (tauR - time + 2.0 * tauS)) -
           tauR * (tauR - time + 2.0 * tauS) * std::asin(std::sqrt((time - 2.0 * tauS) / tauR)) -
           0.75 * tauR * tauR *
               std::atan(std::sqrt((tauR - time + 2.0 * tauS) / (time - 2.0 * tauS)));
  };

  auto c5 = [&]() { return 0.5 * M_PI * tauR * (time - tauR); };

  auto c6 = [&]() { return 0.5 * M_PI * tauR * (2.0 * tauS - time + tauR); };

  if (tauR > 2.0 * tauS) {
    if (time <= 0) {
      return 0;
    } else if (time <= tauS) {
      return k * (c1() + c2());
    } else if (time <= 2.0 * tauS) {
      return k * (c1() - c2() + c3());
    } else if (time < tauR) {
      return k * (c1() + c3() + c4());
    } else if (time < tauR + tauS) {
      return k * (c3() + c4() + c5());
    } else if (time < tauR + 2.0 * tauS) {
      return k * (c4() + c6());
    } else {
      return 0;
    }
  } else {
    if (time <= 0) {
      return 0;
    } else if (time <= tauS) {
      return k * (c1() + c2());
    } else if (time < tauR) {
      return k * (c1() - c2() + c3());
    } else if (time <= 2.0 * tauS) {
      return k * (c5() + c3() - c2());
    } else if (time < tauR + tauS) {
      return k * (c3() + c4() + c5());
    } else if (time < tauR + 2.0 * tauS) {
      return k * (c4() + c6());
    } else {
      return 0;
    }
  }
}
} // namespace seissol::regularizedYoffe

#endif // SEISSOL_REGULARIZEDYOFFE_H

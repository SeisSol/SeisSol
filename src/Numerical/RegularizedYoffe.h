// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_
#define SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_

#include "Common/Marker.h"

#include <cmath>

namespace seissol::regularizedYoffe {
/**
 * Implementation of the regularized Yoffe function defined in Appendix of Tinti et al. (2005)
 *
 * The regularization is a convolution with a triangle of half width tauS, and the closed form
 * below is the corresponding second difference: c3 is -2 c1 shifted by tauS and c4 is c1
 * shifted by 2 tauS, with c2, c5 and c6 the corrections where a shifted argument leaves the
 * support. A second difference of step tauS over a function of scale tauR cancels down to a
 * result of relative size (tauS / tauR)^2, so the terms lose (tauR / tauS)^2 in relative
 * accuracy no matter how they are arranged. Measured against the convolution at 40 digits,
 * relative to the peak of the function, for tauR / tauS from 1 to 200: 7e-16 to 1e-12 in
 * double and 3e-07 to 7e-04 in single.
 *
 * Compute is where that sum is carried. It stays T by default, so that single precision stays
 * single precision all the way through, which is what hardware with a crippled double rate
 * needs. Passing double instead pins the error at the final rounding, 3e-08 relative to the
 * peak whatever the ratio, and is worth its cost where tauR / tauS is extreme.
 */
template <typename T, typename Compute = T>
SEISSOL_HOSTDEVICE inline T regularizedYoffe(T timeIn, T tauSIn, T tauRIn) {
  const auto time = static_cast<Compute>(timeIn);
  const auto tauS = static_cast<Compute>(tauSIn);
  const auto tauR = static_cast<Compute>(tauRIn);
  // constants to respect the target precision
  constexpr Compute K2 = 2.0;
  constexpr Compute K025 = 0.25;
  constexpr Compute K0375 = 0.375;
  constexpr Compute K05 = 0.5;
  constexpr Compute K075 = 0.75;
  constexpr Compute K15 = 1.5;
  constexpr Compute KPi = M_PI;

  const auto k = K2 / (KPi * tauR * tauS * tauS);
  // c1 to c6 are analytical functions used for building the regularized Yoffe function
  const auto c1 = [&]() {
    return (K05 * time + K025 * tauR) * std::sqrt(time * (tauR - time)) +
           (time * tauR - tauR * tauR) * std::asin(std::sqrt(time / tauR)) -
           K075 * tauR * tauR * std::atan(std::sqrt((tauR - time) / time));
  };

  const auto c2 = [&] { return K0375 * KPi * tauR * tauR; };

  const auto c3 = [&]() {
    return (tauS - time - K05 * tauR) * std::sqrt((time - tauS) * (tauR - time + tauS)) +
           tauR * (K2 * tauR - K2 * time + K2 * tauS) * std::asin(std::sqrt((time - tauS) / tauR)) +
           K15 * tauR * tauR * std::atan(std::sqrt((tauR - time + tauS) / (time - tauS)));
  };

  const auto c4 = [&]() {
    // 2 typos fixed in the second term compared with Tinti et al. 2005
    return (-tauS + K05 * time + K025 * tauR) *
               std::sqrt((time - K2 * tauS) * (tauR - time + K2 * tauS)) -
           tauR * (tauR - time + K2 * tauS) * std::asin(std::sqrt((time - K2 * tauS) / tauR)) -
           K075 * tauR * tauR *
               std::atan(std::sqrt((tauR - time + K2 * tauS) / (time - K2 * tauS)));
  };

  const auto c5 = [&]() { return K05 * KPi * tauR * (time - tauR); };

  const auto c6 = [&]() { return K05 * KPi * tauR * (K2 * tauS - time + tauR); };

  if (tauR > K2 * tauS) {
    if (time <= 0) {
      return static_cast<T>(0);
    } else if (time <= tauS) {
      return static_cast<T>(k * (c1() + c2()));
    } else if (time <= K2 * tauS) {
      return static_cast<T>(k * (c1() - c2() + c3()));
    } else if (time < tauR) {
      return static_cast<T>(k * (c1() + c3() + c4()));
    } else if (time < tauR + tauS) {
      return static_cast<T>(k * (c3() + c4() + c5()));
    } else if (time < tauR + K2 * tauS) {
      return static_cast<T>(k * (c4() + c6()));
    } else {
      return static_cast<T>(0);
    }
  } else {
    if (time <= 0) {
      return static_cast<T>(0);
    } else if (time <= tauS) {
      return static_cast<T>(k * (c1() + c2()));
    } else if (time < tauR) {
      return static_cast<T>(k * (c1() - c2() + c3()));
    } else if (time <= K2 * tauS) {
      return static_cast<T>(k * (c5() + c3() - c2()));
    } else if (time < tauR + tauS) {
      return static_cast<T>(k * (c3() + c4() + c5()));
    } else if (time < tauR + K2 * tauS) {
      return static_cast<T>(k * (c4() + c6()));
    } else {
      return static_cast<T>(0);
    }
  }
}
} // namespace seissol::regularizedYoffe

#endif // SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_

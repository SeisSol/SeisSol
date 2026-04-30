// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_
#define SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_

#include "Common/Marker.h"

namespace seissol::regularizedYoffe {
/**
 * Implementation of the regularized Yoffe function defined in Appendix of Tinti et al. (2005)
 */
template <typename T>
SEISSOL_HOSTDEVICE inline T regularizedYoffe(T time, T tauS, T tauR) {
  // constants to respect the target precision
  constexpr T K2 = 2.0;
  constexpr T K025 = 0.25;
  constexpr T K0375 = 0.375;
  constexpr T K05 = 0.5;
  constexpr T K075 = 0.75;
  constexpr T K15 = 1.5;
  constexpr T KPi = M_PI;

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
      return 0;
    } else if (time <= tauS) {
      return k * (c1() + c2());
    } else if (time <= K2 * tauS) {
      return k * (c1() - c2() + c3());
    } else if (time < tauR) {
      return k * (c1() + c3() + c4());
    } else if (time < tauR + tauS) {
      return k * (c3() + c4() + c5());
    } else if (time < tauR + K2 * tauS) {
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
    } else if (time <= K2 * tauS) {
      return k * (c5() + c3() - c2());
    } else if (time < tauR + tauS) {
      return k * (c3() + c4() + c5());
    } else if (time < tauR + K2 * tauS) {
      return k * (c4() + c6());
    } else {
      return 0;
    }
  }
}
} // namespace seissol::regularizedYoffe

#endif // SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_

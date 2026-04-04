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
  constexpr T k2 = 2.0;
  constexpr T k025 = 0.25;
  constexpr T k0375 = 0.375;
  constexpr T k05 = 0.5;
  constexpr T k075 = 0.75;
  constexpr T k15 = 1.5;
  constexpr T kPi = M_PI;

  const auto k = k2 / (kPi * tauR * tauS * tauS);
  // c1 to c6 are analytical functions used for building the regularized Yoffe function
  const auto c1 = [&]() {
    return (k05 * time + k025 * tauR) * std::sqrt(time * (tauR - time)) +
           (time * tauR - tauR * tauR) * std::asin(std::sqrt(time / tauR)) -
           k075 * tauR * tauR * std::atan(std::sqrt((tauR - time) / time));
  };

  const auto c2 = [&] { return k0375 * kPi * tauR * tauR; };

  const auto c3 = [&]() {
    return (tauS - time - k05 * tauR) * std::sqrt((time - tauS) * (tauR - time + tauS)) +
           tauR * (k2 * tauR - k2 * time + k2 * tauS) * std::asin(std::sqrt((time - tauS) / tauR)) +
           k15 * tauR * tauR * std::atan(std::sqrt((tauR - time + tauS) / (time - tauS)));
  };

  const auto c4 = [&]() {
    // 2 typos fixed in the second term compared with Tinti et al. 2005
    return (-tauS + k05 * time + k025 * tauR) *
               std::sqrt((time - k2 * tauS) * (tauR - time + k2 * tauS)) -
           tauR * (tauR - time + k2 * tauS) * std::asin(std::sqrt((time - k2 * tauS) / tauR)) -
           k075 * tauR * tauR *
               std::atan(std::sqrt((tauR - time + k2 * tauS) / (time - k2 * tauS)));
  };

  const auto c5 = [&]() { return k05 * kPi * tauR * (time - tauR); };

  const auto c6 = [&]() { return k05 * kPi * tauR * (k2 * tauS - time + tauR); };

  if (tauR > k2 * tauS) {
    if (time <= 0) {
      return 0;
    } else if (time <= tauS) {
      return k * (c1() + c2());
    } else if (time <= k2 * tauS) {
      return k * (c1() - c2() + c3());
    } else if (time < tauR) {
      return k * (c1() + c3() + c4());
    } else if (time < tauR + tauS) {
      return k * (c3() + c4() + c5());
    } else if (time < tauR + k2 * tauS) {
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
    } else if (time <= k2 * tauS) {
      return k * (c5() + c3() - c2());
    } else if (time < tauR + tauS) {
      return k * (c3() + c4() + c5());
    } else if (time < tauR + k2 * tauS) {
      return k * (c4() + c6());
    } else {
      return 0;
    }
  }
}
} // namespace seissol::regularizedYoffe

#endif // SEISSOL_SRC_NUMERICAL_REGULARIZEDYOFFE_H_

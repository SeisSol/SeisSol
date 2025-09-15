// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_

#include <cstdint>
#include <limits>
#include <type_traits>

namespace seissol::dr::friction_law::rs {
// If the SR is too close to zero, we will have problems (NaN)
// as a consequence, the SR is affected the AlmostZero value when too small
// For double precision 1e-45 is a chosen by trial and error. For single precision, this value is
// too small, so we use 1e-35
template <typename RealT>
constexpr RealT almostZero() {
  if constexpr (std::is_same<RealT, double>()) {
    return 1e-45;
  } else if constexpr (std::is_same<RealT, float>()) {
    return 1e-35;
  } else {
    return std::numeric_limits<RealT>::min();
  }
}

struct Settings {
  /**
   * Parameters of the optimisation loops
   * absolute tolerance on the function to be optimized
   * This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be
   * the most adapted Number of iteration in the loops
   */

  const uint32_t maxNumberSlipRateUpdates{60};
  const uint32_t numberStateVariableUpdates{2};
  const double newtonTolerance{1e-8};
};
} // namespace seissol::dr::friction_law::rs

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_RATEANDSTATECOMMON_H_

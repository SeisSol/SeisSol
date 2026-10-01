// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_TPCOMMON_T_H_
#define SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_TPCOMMON_T_H_

#include <doctest.h>

#include "DynamicRupture/FrictionLaws/TPCommon.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "TestHelper.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace seissol::unit_test {

namespace tptest {

namespace tp = seissol::dr::friction_law::tp;
namespace misc = seissol::dr::misc;

/// The grid is a runtime parameter now, so every property has to hold at whatever count the
/// parameter file names -- not only at the default.
inline const std::vector<std::size_t>& counts() {
  static const std::vector<std::size_t> Values = {
      2, 3, 8, 17, 32, misc::DefaultTpGridPoints, 100, 257};
  return Values;
}

/**
 * The peak of the standard normal density. The inverse transform of the Gaussian shear-zone
 * spectrum has to reproduce exactly this at z = 0.
 */
const double GaussianPeak = 1.0 / std::sqrt(2.0 * M_PI);

} // namespace tptest

// ---------------------------------------------------------------------------
// GridPoints
// ---------------------------------------------------------------------------

TEST_CASE("DR TPCommon GridPoints" * doctest::test_suite("dynamicrupture")) {
  using namespace tptest;

  SUBCASE("As long as the count says") {
    for (const std::size_t count : counts()) {
      CAPTURE(count);
      const tp::GridPoints<double> gridPoints(count);
      CHECK(gridPoints.size() == count);
      CHECK(gridPoints.data().size() == count);
    }
  }

  SUBCASE("Strictly increasing and positive") {
    for (const std::size_t count : counts()) {
      const tp::GridPoints<double> gridPoints(count);
      for (std::size_t i = 0; i < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        CHECK(gridPoints[i] > 0.0);
        if (i > 0) {
          CHECK(gridPoints[i] > gridPoints[i - 1]);
        }
      }
    }
  }

  SUBCASE("A geometric progression of ratio exp(TpLogDz)") {
    // the logarithmic grid of Noda & Lapusta (2010), eq. (14). The quadrature weights in
    // InverseFourierCoefficients assume exactly this spacing.
    const double ratio = std::exp(misc::TpLogDz);
    for (const std::size_t count : counts()) {
      const tp::GridPoints<double> gridPoints(count);
      for (std::size_t i = 1; i < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        CHECK(gridPoints[i] / gridPoints[i - 1] == doctest::Approx(ratio).epsilon(1e-14));
      }
    }
  }

  SUBCASE("Anchored at TpMaxWaveNumber, whatever the count") {
    // the count is a runtime parameter, so the grid has to end at the same wave number for every
    // one of them and simply resolve it more or less finely. A count baked into the exponent
    // instead would leave the largest wave number of a 32-point grid at 2e-3 and take every
    // derived quantity with it.
    for (const std::size_t count : counts()) {
      CAPTURE(count);
      const tp::GridPoints<double> gridPoints(count);
      CHECK(gridPoints[count - 1] == doctest::Approx(misc::TpMaxWaveNumber));
      CHECK(gridPoints[0] ==
            doctest::Approx(misc::TpMaxWaveNumber *
                            std::exp(-misc::TpLogDz * static_cast<double>(count - 1))));
    }
  }

  SUBCASE("data() exposes the same values") {
    for (const std::size_t count : counts()) {
      const tp::GridPoints<double> gridPoints(count);
      const auto& values = gridPoints.data();
      for (std::size_t i = 0; i < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        CHECK(values[i] == gridPoints[i]);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// InverseFourierCoefficients
// ---------------------------------------------------------------------------

TEST_CASE("DR TPCommon InverseFourierCoefficients" * doctest::test_suite("dynamicrupture")) {
  using namespace tptest;
  const double norm = std::sqrt(2.0 / M_PI);

  SUBCASE("As long as the count says") {
    for (const std::size_t count : counts()) {
      CAPTURE(count);
      CHECK(tp::InverseFourierCoefficients<double>(count).size() == count);
    }
  }

  SUBCASE("Interior weights are sqrt(2/pi) l_i TpLogDz") {
    // dl = l dlog(l) on a logarithmic grid, so the trapezoidal weight carries a factor l
    for (const std::size_t count : counts()) {
      const tp::GridPoints<double> gridPoints(count);
      const tp::InverseFourierCoefficients<double> coefficients(count);
      for (std::size_t i = 1; i + 1 < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        CHECK(coefficients[i] == doctest::Approx(norm * gridPoints[i] * misc::TpLogDz));
      }
    }
  }

  SUBCASE("Boundary weights are special-cased") {
    for (const std::size_t count : counts()) {
      CAPTURE(count);
      const tp::GridPoints<double> gridPoints(count);
      const tp::InverseFourierCoefficients<double> coefficients(count);
      CHECK(coefficients[0] == doctest::Approx(norm * gridPoints[0] * (1.0 + misc::TpLogDz)));
      if (count > 1) {
        CHECK(coefficients[count - 1] ==
              doctest::Approx(norm * gridPoints[count - 1] * 0.5 * misc::TpLogDz));
      }
    }
  }

  SUBCASE("All weights positive") {
    for (const std::size_t count : counts()) {
      const tp::InverseFourierCoefficients<double> coefficients(count);
      for (std::size_t i = 0; i < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        CHECK(coefficients[i] > 0.0);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// GaussianHeatSource
// ---------------------------------------------------------------------------

TEST_CASE("DR TPCommon GaussianHeatSource" * doctest::test_suite("dynamicrupture")) {
  using namespace tptest;

  SUBCASE("As long as the count says") {
    for (const std::size_t count : counts()) {
      CAPTURE(count);
      CHECK(tp::GaussianHeatSource<double>(count).size() == count);
    }
  }

  SUBCASE("Equals exp(-l^2/2) / sqrt(2 pi)") {
    for (const std::size_t count : counts()) {
      const tp::GridPoints<double> gridPoints(count);
      const tp::GaussianHeatSource<double> heatSource(count);
      for (std::size_t i = 0; i < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        const double expected =
            std::exp(-0.5 * gridPoints[i] * gridPoints[i]) / std::sqrt(2.0 * M_PI);
        CHECK(heatSource[i] == doctest::Approx(expected));
      }
    }
  }

  SUBCASE("Decreasing, bounded by the Gaussian peak") {
    for (const std::size_t count : counts()) {
      const tp::GaussianHeatSource<double> heatSource(count);
      for (std::size_t i = 0; i < count; ++i) {
        CAPTURE(count);
        CAPTURE(i);
        CHECK(heatSource[i] > 0.0);
        CHECK(heatSource[i] <= GaussianPeak);
        if (i > 0) {
          CHECK(heatSource[i] <= heatSource[i - 1]);
        }
      }
    }
  }

  SUBCASE("The far end has decayed by exp(-l^2/2) at TpMaxWaveNumber") {
    for (const std::size_t count : counts()) {
      CAPTURE(count);
      const tp::GaussianHeatSource<double> heatSource(count);
      CHECK(heatSource[count - 1] ==
            doctest::Approx(std::exp(-0.5 * misc::TpMaxWaveNumber * misc::TpMaxWaveNumber) *
                            GaussianPeak));
    }
  }

  SUBCASE("A fine grid starts at the peak") {
    // the smallest wave number of the default grid is 2e-7, so the spectrum is still flat there
    const tp::GaussianHeatSource<double> heatSource(misc::DefaultTpGridPoints);
    CHECK(heatSource[0] == doctest::Approx(GaussianPeak));
  }
}

// ---------------------------------------------------------------------------
// the discretisation as a whole
// ---------------------------------------------------------------------------

TEST_CASE("DR TPCommon inverse transform reproduces the Gaussian peak" *
          doctest::test_suite("dynamicrupture")) {
  using namespace tptest;
  // The weights implement the inverse cosine transform of Noda & Lapusta (2010),
  //     Theta(z) = sqrt(2/pi) int_0^inf Theta_hat(l) cos(l z) dl,
  // discretised on the logarithmic grid, dl = l TpLogDz. With the Gaussian shear-zone spectrum
  // Theta_hat(l) = exp(-l^2/2) / sqrt(2 pi) this evaluates at z = 0 to
  //     sum_i c_i h_i = 1 / sqrt(2 pi) = 0.39894228...
  //
  // This one number ties the grid, the quadrature weights and the heat source together, and it is
  // the sharpest check available on the discretisation: changing TpLogDz or TpMaxWaveNumber without
  // re-deriving the weights breaks it at once. It is also what says how many grid points a run
  // actually needs, which is now a parameter.
  const auto quadrature = [](std::size_t count) {
    const tp::InverseFourierCoefficients<double> coefficients(count);
    const tp::GaussianHeatSource<double> heatSource(count);
    double sum = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
      sum += coefficients[i] * heatSource[i];
    }
    return sum;
  };

  SUBCASE("At the default count") {
    // measured with the current constants: 4.1e-9 relative
    CHECK(quadrature(misc::DefaultTpGridPoints) == doctest::Approx(GaussianPeak).epsilon(1e-7));
  }

  SUBCASE("Converging as the count grows") {
    // the grid ends at TpMaxWaveNumber for every count, so a larger count refines the spacing
    // rather than extending the range, and the quadrature has to improve monotonically far enough
    // to make a parameter worth having
    const double coarse = std::abs(quadrature(20) / GaussianPeak - 1.0);
    const double medium = std::abs(quadrature(40) / GaussianPeak - 1.0);
    const double fine = std::abs(quadrature(misc::DefaultTpGridPoints) / GaussianPeak - 1.0);
    CAPTURE(coarse);
    CAPTURE(medium);
    CAPTURE(fine);
    CHECK(medium < coarse);
    CHECK(fine < medium);
  }

  SUBCASE("A count too small to resolve the spectrum is visibly wrong") {
    // which is the other half of making it a parameter: the error has to be there to be seen
    CHECK(std::abs(quadrature(8) / GaussianPeak - 1.0) > 1e-3);
  }
}

TEST_CASE("DR TPCommon in working precision" * doctest::test_suite("dynamicrupture")) {
  using namespace tptest;
  // The tables are built in double and stored as `real`. The pressurization divides the
  // coefficients by the half width of the shear zone and multiplies the heat source by tau V, so
  // neither table may arrive as a subnormal.
  const auto count = static_cast<std::size_t>(misc::DefaultTpGridPoints);
  const tp::GridPoints<> gridPoints(count);
  const tp::InverseFourierCoefficients<> coefficients(count);
  const tp::GaussianHeatSource<> heatSource(count);

  SUBCASE("No underflow or overflow in the stored tables") {
    for (std::size_t i = 0; i < count; ++i) {
      CAPTURE(i);
      CHECK(std::isnormal(gridPoints[i]));
      CHECK(std::isnormal(coefficients[i]));
      CHECK(std::isnormal(heatSource[i]));
    }
  }

  SUBCASE("The quadrature identity survives the narrowing") {
    double sum = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
      sum += static_cast<double>(coefficients[i]) * static_cast<double>(heatSource[i]);
    }
    // measured: 4.1e-9 relative for a double build, 3.9e-9 for a single-precision one
    CHECK(sum == doctest::Approx(GaussianPeak).epsilon(1e-5));
  }
}

} // namespace seissol::unit_test

#endif // SEISSOL_TESTS_DYNAMICRUPTURE_FRICTIONLAWS_TPCOMMON_T_H_

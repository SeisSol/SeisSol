// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Kernels/PointSourceCluster.h"
#include "Kernels/Precision.h"

#include "doctest.h"

#include <cmath>
#include <vector>

namespace seissol::unit_test {
TEST_CASE("computeSampleTimeIntegral") {
  constexpr double Pi = 3.14159265358979323846264338327950;

  std::size_t sampleSize = 1000;
  double samplingInterval = 1.0 / (sampleSize - 1);
  auto sr = std::vector<real>(sampleSize);
  for (std::size_t i = 0; i < sampleSize; ++i) {
    sr[i] = std::sin(Pi * i * samplingInterval);
  }

  const auto computeIntegral = [&](double from, double to, double onsetTime) {
    return kernels::computeSampleTimeIntegral(
        from, to, onsetTime, samplingInterval, sr.data(), sampleSize);
  };

  const auto ref = [&](double from, double to) {
    return (-std::cos(Pi * to) + std::cos(Pi * from)) / Pi;
  };

  // Basic case
  CHECK(computeIntegral(0.0, 1.0, 0.0) == doctest::Approx(ref(0.0, 1.0)));
  // Onset time
  CHECK(computeIntegral(42.1005, 42.8995, 42.0) == doctest::Approx(ref(0.1005, 0.8995)));

  // Partial or no overlap
  CHECK(computeIntegral(-1.0, 0.5, 0.0) == doctest::Approx(ref(0.0, 0.5)));
  CHECK(computeIntegral(0.75, 1.5, 0.0) == doctest::Approx(ref(0.75, 1.0)));
  CHECK(computeIntegral(1.0, 3.0, 0.0) == 0.0);

  // Involves only two samples
  CHECK(computeIntegral(2.1 * samplingInterval, 2.8 * samplingInterval, 0.0) ==
        doctest::Approx(ref(2.1 * samplingInterval, 2.8 * samplingInterval)));
}

} // namespace seissol::unit_test

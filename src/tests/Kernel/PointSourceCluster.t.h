// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include <Kernels/PointSourceCluster.h>
#include <Kernels/precision.hpp>

#include "doctest.h"

#include <vector>
#include <cmath>

namespace seissol::unit_test {
TEST_CASE("computeSampleTimeIntegral") {
  constexpr double pi = 3.14159265358979323846264338327950;

  std::size_t sampleSize = 1000;
  double samplingInterval = 1.0 / (sampleSize - 1);
  auto sr = std::vector<real>(sampleSize);
  for (std::size_t i = 0; i < sampleSize; ++i) {
    sr[i] = std::sin(pi * i * samplingInterval);
  }

  auto const computeIntegral = [&](double from, double to, double onsetTime) {
    return kernels::computeSampleTimeIntegral(
        from, to, onsetTime, samplingInterval, sr.data(), sampleSize);
  };

  auto const ref = [&](double from, double to) {
    return (-std::cos(pi * to) + std::cos(pi * from)) / pi;
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

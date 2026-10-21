// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/DeviceTimeline.h"

#include <utility>
#include <vector>

namespace seissol::unit_test {

TEST_CASE("The busy time of the device counts each moment once" *
          doctest::test_suite("monitoring")) {
  using Intervals = std::vector<std::pair<double, double>>;
  CHECK(DeviceTimeline::unionLength({}) == doctest::Approx(0.0));
  CHECK(DeviceTimeline::unionLength(Intervals{{1.0, 3.0}}) == doctest::Approx(2.0));
  // apart, overlapping, contained, and touching; in any order
  CHECK(DeviceTimeline::unionLength(Intervals{{5.0, 6.0}, {1.0, 3.0}}) == doctest::Approx(3.0));
  CHECK(DeviceTimeline::unionLength(Intervals{{2.0, 4.0}, {1.0, 3.0}}) == doctest::Approx(3.0));
  CHECK(DeviceTimeline::unionLength(Intervals{{1.0, 5.0}, {2.0, 3.0}}) == doctest::Approx(4.0));
  CHECK(DeviceTimeline::unionLength(Intervals{{3.0, 4.0}, {1.0, 3.0}}) == doctest::Approx(3.0));
  CHECK(DeviceTimeline::unionLength(Intervals{{0.0, 1.0}, {0.5, 2.0}, {1.5, 3.0}, {4.0, 5.0}}) ==
        doctest::Approx(4.0));
}

TEST_CASE("Without timing on the device, an epoch has no device time" *
          doctest::test_suite("monitoring")) {
  if (DeviceTimeline::active()) {
    return;
  }
  DeviceTimeline timeline;
  timeline.beginEpoch();
  const auto activity = timeline.endEpoch();
  CHECK(activity.span == doctest::Approx(0.0));
  CHECK(activity.busy == doctest::Approx(0.0));
  CHECK(activity.computing == doctest::Approx(0.0));
}

} // namespace seissol::unit_test

// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <Common/Iterator.h>

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Pre-C++20/23 Iterators") {
  std::vector<int32_t> v1{2, 4, 6, 8, 10};
  std::vector<float> v2{.3, .6, .9, .12, .15};

  std::vector<std::tuple<int32_t, float>> v12;
  v12.reserve(v1.size());
  for (std::size_t i = 0; i < v1.size(); ++i) {
    v12.emplace_back(v1[i], v2[i]);
  }

  std::size_t n4 = 0;
  for (std::size_t i = 0; i < v1.size(); ++i) {
    if (v1[i] % 4 == 0) {
      ++n4;
    }
  }

  SUBCASE("Zip") {
    std::size_t i = 0;
    for (const auto t : common::zip(v1, v2)) {
      REQUIRE(t == v12[i]);
      ++i;
    }

    REQUIRE(i == v12.size());
  }

  SUBCASE("Enumerate") {
    std::size_t j = 0;
    for (const auto [i, e] : common::enumerate(v1)) {
      REQUIRE(i == j);
      REQUIRE(e == v1[j]);
      ++j;
    }
  }

  SUBCASE("Filter") {
    // no short-cut method here (yet, or ever---it'll be gone once we switch to C++20 anyways)

    auto it =
        common::FilteredIterator(v1.begin(), v1.end(), [](auto value) { return value % 4 == 0; });
    auto itEnd =
        common::FilteredIterator(v1.end(), v1.end(), [](auto value) { return value % 4 == 0; });
    std::size_t i = 0;
    for (; it != itEnd; ++it) {
      REQUIRE(*it % 4 == 0);
      ++i;
    }

    REQUIRE(i == n4);
  }
}

} // namespace seissol::unit_test

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Dimension.h"

#include <cstdint>
#include <vector>

namespace seissol::unit_test {

namespace dimensiontest {
using namespace seissol::io::writer;
} // namespace dimensiontest

using namespace dimensiontest;

TEST_CASE("IO/Dimension: a shape says how it joins the dataset" * doctest::test_suite("io")) {
  SUBCASE("The common shape has a leading distributed dimension and nothing that grows") {
    const auto dimensions = makeDimensions({3}, true);
    REQUIRE(dimensions.size() == 2);
    CHECK(dimensions[0].isDistributed());
    CHECK_FALSE(dimensions[0].isAppended());
    CHECK_FALSE(dimensions[1].isDistributed());
    CHECK(dimensions[1].size == 3);
  }

  SUBCASE("An appended dimension carries what one write contributes") {
    const auto dimension = Dimension::appended(16);
    CHECK(dimension.isAppended());
    CHECK_FALSE(dimension.isDistributed());
    CHECK(dimension.size == 16);
  }

  SUBCASE("A flat array is distributed and appended at once") {
    const auto dimension = Dimension::distributedAppended();
    CHECK(dimension.isAppended());
    CHECK(dimension.isDistributed());
  }
}

TEST_CASE("IO/Dimension: a source keeps the shape it was given" * doctest::test_suite("io")) {
  const std::vector<std::int64_t> values{1, 2, 3, 4, 5, 6};
  // what a receiver table looks like: samples grow, receivers are split across the ranks, and
  // every entry holds the same quantities
  const std::vector<Dimension> dimensions{
      Dimension::appended(2), Dimension::distributed(), Dimension::replicated(3)};
  const auto source = WriteInline::createShaped<std::int64_t>(dimensions, values);

  REQUIRE(source->dimensions().size() == 3);
  CHECK(source->dimensions()[0].isAppended());
  CHECK(source->dimensions()[1].isDistributed());
  CHECK(source->distributed());
  // only what an entry holds, so neither the growing nor the distributed dimension
  CHECK(source->shape() == std::vector<std::size_t>{3});

  SUBCASE("and survives serialization") {
    const auto restored = DataSource::deserialize(source->serialize());
    REQUIRE(restored->dimensions().size() == 3);
    CHECK(restored->dimensions()[0].isAppended());
    CHECK(restored->dimensions()[0].size == 2);
    CHECK(restored->dimensions()[1].isDistributed());
    CHECK(restored->dimensions()[2].size == 3);
    CHECK(restored->shape() == std::vector<std::size_t>{3});
  }
}

} // namespace seissol::unit_test

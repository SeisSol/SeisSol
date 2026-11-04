// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/RunId.h"

#include <set>
#include <string>

namespace seissol::unit_test {

TEST_CASE("Run IDs" * doctest::test_suite("common")) {
  SUBCASE("a new run ID is a version 4 UUID") {
    const auto runId = seissol::newRunId();
    REQUIRE(runId.size() == seissol::RunIdLength);
    CHECK(runId[8] == '-');
    CHECK(runId[13] == '-');
    CHECK(runId[18] == '-');
    CHECK(runId[23] == '-');
    // version and variant
    CHECK(runId[14] == '4');
    CHECK(std::string("89ab").find(runId[19]) != std::string::npos);
    CHECK(seissol::parseRunId(runId) == runId);
  }

  SUBCASE("new run IDs differ") {
    std::set<std::string> runIds;
    for (int i = 0; i < 100; ++i) {
      runIds.insert(seissol::newRunId());
    }
    CHECK(runIds.size() == 100);
  }

  SUBCASE("parsing normalizes to lower case") {
    CHECK(seissol::parseRunId("123E4567-E89B-42D3-A456-426614174000") ==
          std::string("123e4567-e89b-42d3-a456-426614174000"));
  }

  SUBCASE("parsing rejects what is not a UUID") {
    CHECK_FALSE(seissol::parseRunId("").has_value());
    CHECK_FALSE(seissol::parseRunId("123e4567e89b42d3a456426614174000").has_value());
    CHECK_FALSE(seissol::parseRunId("123e4567-e89b-42d3-a456-42661417400g").has_value());
    CHECK_FALSE(seissol::parseRunId("123e4567-e89b-42d3-a456_426614174000").has_value());
    CHECK_FALSE(seissol::parseRunId("123e4567-e89b-42d3-a456-4266141740001").has_value());
  }
}

} // namespace seissol::unit_test

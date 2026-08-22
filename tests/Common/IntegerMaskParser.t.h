// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/IntegerMaskParser.h"

namespace seissol::unit_test {

using namespace seissol;

TEST_CASE("Integer Mask Parser" * doctest::test_suite("common")) {
  SUBCASE("Parse a list of integers") {
    auto parsedResult = IntegerMaskParser::parse("10,20,30");
    CHECK(parsedResult.has_value());
    auto result = *parsedResult;
    CHECK(result.size() == 3);
    CHECK(result[0][0] == 10);
    CHECK(result[1][0] == 20);
    CHECK(result[2][0] == 30);
  }

  SUBCASE("Parse a list of ranges") {
    auto parsedResult = IntegerMaskParser::parse("0-2,3-5,6-8");
    CHECK(parsedResult.has_value());
    auto result = *parsedResult;
    CHECK(result.size() == 3);
    CHECK(result[0] == std::vector<int>{0, 1, 2});
    CHECK(result[1] == std::vector<int>{3, 4, 5});
    CHECK(result[2] == std::vector<int>{6, 7, 8});
  }

  SUBCASE("Parse a mixture of allowed formats") {
    auto parsedResult = IntegerMaskParser::parse("0,{3,4,5},6-8");
    CHECK(parsedResult.has_value());
    auto result = *parsedResult;
    CHECK(result.size() == 3);
    CHECK(result[0] == std::vector<int>{0});
    CHECK(result[1] == std::vector<int>{3, 4, 5});
    CHECK(result[2] == std::vector<int>{6, 7, 8});
  }

  SUBCASE("Parse with the leading error / integer") {
    auto parsedResult = IntegerMaskParser::parse("0th,{3,4,5},6-8");
    CHECK(!parsedResult.has_value());
  }

  SUBCASE("Parse with an error at the middle / list") {
    auto parsedResult = IntegerMaskParser::parse("0-2,{3,x4,5},6-8");
    CHECK(parsedResult.has_value());
    auto result = *parsedResult;
    CHECK(result.size() == 1);
    CHECK(result[0] == std::vector<int>{0, 1, 2});
  }

  SUBCASE("Parse with an error at the end / range") {
    auto parsedResult = IntegerMaskParser::parse("0-2,{3,4,5},6--8");
    CHECK(parsedResult.has_value());
    auto result = *parsedResult;
    CHECK(result.size() == 2);
    CHECK(result[0] == std::vector<int>{0, 1, 2});
    CHECK(result[1] == std::vector<int>{3, 4, 5});
  }
}
} // namespace seissol::unit_test

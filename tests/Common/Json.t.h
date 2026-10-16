// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Filesystem.h"
#include "Common/Json.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::unit_test {

using seissol::json::Value;

TEST_CASE("JSON scalars" * doctest::test_suite("common")) {
  CHECK(Value().dump() == "null");
  CHECK(Value(nullptr).dump() == "null");
  CHECK(Value(true).dump() == "true");
  CHECK(Value(false).dump() == "false");
  CHECK(Value(42).dump() == "42");
  CHECK(Value(-7).dump() == "-7");
  CHECK(Value(std::numeric_limits<std::uint64_t>::max()).dump() == "18446744073709551615");
  CHECK(Value(std::numeric_limits<std::int64_t>::min()).dump() == "-9223372036854775808");
  CHECK(Value("text").dump() == "\"text\"");
  CHECK(Value(std::string("text")).dump() == "\"text\"");
}

TEST_CASE("JSON numbers" * doctest::test_suite("common")) {
  SUBCASE("shortest spelling that reads back exactly") {
    CHECK(Value(0.1).dump() == "0.1");
    CHECK(Value(1.5).dump() == "1.5");
    CHECK(Value(1e-300).dump() == "1e-300");
    for (const double value :
         {1.0 / 3.0, 2.0 / 3.0, 123456.789e12, 5e-324, 1.7976931348623157e308}) {
      const auto text = Value(value).dump();
      CHECK(std::strtod(text.c_str(), nullptr) == value);
    }
  }
  SUBCASE("not finite numbers are null") {
    CHECK(Value(std::numeric_limits<double>::infinity()).dump() == "null");
    CHECK(Value(-std::numeric_limits<double>::infinity()).dump() == "null");
    CHECK(Value(std::numeric_limits<double>::quiet_NaN()).dump() == "null");
  }
}

TEST_CASE("JSON string escaping" * doctest::test_suite("common")) {
  CHECK(json::quote("a\"b") == "\"a\\\"b\"");
  CHECK(json::quote("a\\b") == "\"a\\\\b\"");
  CHECK(json::quote("line\nnext\ttab\r") == "\"line\\nnext\\ttab\\r\"");
  CHECK(json::quote(std::string("\x01\x1f", 2)) == "\"\\u0001\\u001f\"");
  // UTF-8 stays as it is
  CHECK(json::quote("µs") == "\"µs\"");
}

TEST_CASE("JSON objects and arrays" * doctest::test_suite("common")) {
  SUBCASE("objects keep the order members were added in") {
    Value value;
    value["zeta"] = 1;
    value["alpha"] = 2;
    value["mid"]["inner"] = "x";
    CHECK(value.dump(-1) == "{\"zeta\":1,\"alpha\":2,\"mid\":{\"inner\":\"x\"}}");
    CHECK(value.size() == 3);
  }

  SUBCASE("assigning a member again replaces it") {
    Value value;
    value["a"] = 1;
    value["a"] = 2;
    CHECK(value.dump(-1) == "{\"a\":2}");
  }

  SUBCASE("find") {
    Value value;
    value["a"] = 1.5;
    REQUIRE(value.find("a") != nullptr);
    CHECK(value.find("a")->asDouble() == 1.5);
    CHECK(value.find("b") == nullptr);
    CHECK(Value(3).find("a") == nullptr);
  }

  SUBCASE("arrays") {
    Value value;
    value.push(1);
    value.push("two");
    value.push(Value::object());
    CHECK(value.dump(-1) == "[1,\"two\",{}]");
    value.at(2)["three"] = 3;
    CHECK(value.at(2).find("three")->asDouble() == 3);
    CHECK(Value::arrayOf(std::vector<double>{1.0, 2.5}).dump(-1) == "[1,2.5]");
  }

  SUBCASE("empty containers") {
    CHECK(Value::array().dump() == "[]");
    CHECK(Value::object().dump() == "{}");
  }

  SUBCASE("indentation") {
    Value value;
    value["a"] = 1;
    value["b"].push(true);
    CHECK(value.dump(2) == "{\n  \"a\": 1,\n  \"b\": [\n    true\n  ]\n}");
  }
}

TEST_CASE("Atomic file writes" * doctest::test_suite("common")) {
  const std::string path = "json-atomic-test.json";
  REQUIRE(seissol::writeFileAtomically(path, "{\"a\": 1}\n"));
  REQUIRE(seissol::writeFileAtomically(path, "{\"a\": 2}\n"));

  const std::ifstream stream(path);
  std::stringstream content;
  content << stream.rdbuf();
  CHECK(content.str() == "{\"a\": 2}\n");
  CHECK_FALSE(seissol::filesystem::exists(path + ".partial"));

  CHECK_FALSE(seissol::writeFileAtomically("this-directory-does-not-exist/file.json", "{}"));

  seissol::filesystem::remove(path);
}

} // namespace seissol::unit_test

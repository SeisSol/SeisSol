// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/Unit.h"
#include "TestHelper.h"

#include <cstdint>
#include <limits>
#include <string>

namespace seissol::unit_test {

TEST_CASE("formatInteger" * doctest::test_suite("monitoring")) {

  SUBCASE("Zero") { CHECK(formatInteger(0) == "0"); }

  SUBCASE("Single digits") { CHECK(formatInteger(7) == "7"); }

  SUBCASE("Below thousand") {
    CHECK(formatInteger(42) == "42");
    CHECK(formatInteger(999) == "999");
  }

  SUBCASE("Thousand boundary") {
    CHECK(formatInteger(1000) == "1'000");
    CHECK(formatInteger(1001) == "1'001");
  }

  SUBCASE("Larger values") {
    CHECK(formatInteger(1000000) == "1'000'000");
    CHECK(formatInteger(1234567) == "1'234'567");
  }

  SUBCASE("Millions with zero-padded groups") { CHECK(formatInteger(1002003) == "1'002'003"); }

  SUBCASE("Max uint64") {
    const auto result = formatInteger(std::numeric_limits<uint64_t>::max());
    // 18'446'744'073'709'551'615
    CHECK(result.substr(0, 2) == "18");
    CHECK(result.size() > 20);
  }
}

TEST_CASE("SIUnit formatPrefix decimal" * doctest::test_suite("monitoring")) {
  const SIUnit unit("FLOP/s", false);

  SUBCASE("Base range") {
    const auto result = unit.formatPrefix(42.0);
    CHECK(result.find("42.") != std::string::npos);
    CHECK(result.find("FLOP/s") != std::string::npos);
  }

  SUBCASE("Kilo") {
    const auto result = unit.formatPrefix(1500.0);
    CHECK(result.find("1.5") != std::string::npos);
    CHECK(result.find("kFLOP/s") != std::string::npos);
  }

  SUBCASE("Mega") {
    const auto result = unit.formatPrefix(2.5e6);
    CHECK(result.find("2.5") != std::string::npos);
    CHECK(result.find("MFLOP/s") != std::string::npos);
  }

  SUBCASE("Giga") {
    const auto result = unit.formatPrefix(3.0e9);
    CHECK(result.find("3.0") != std::string::npos);
    CHECK(result.find("GFLOP/s") != std::string::npos);
  }

  SUBCASE("Milli") {
    const auto result = unit.formatPrefix(0.005);
    CHECK(result.find("5.0") != std::string::npos);
    CHECK(result.find("mFLOP/s") != std::string::npos);
  }

  SUBCASE("Negative value") {
    const auto result = unit.formatPrefix(-1500.0);
    CHECK(result.find("-1.5") != std::string::npos);
    CHECK(result.find("kFLOP/s") != std::string::npos);
  }

  SUBCASE("Zero") {
    const auto result = unit.formatPrefix(0.0);
    CHECK(result.find("0.0") != std::string::npos);
    CHECK(result.find("FLOP/s") != std::string::npos);
  }

  SUBCASE("With error bar") {
    const auto result = unit.formatPrefix(1500.0, 100.0);
    CHECK(result.find('(') != std::string::npos);
    CHECK(result.find(')') != std::string::npos);
  }
}

TEST_CASE("SIUnit formatPrefix binary" * doctest::test_suite("monitoring")) {
  const SIUnit unit("B", true);

  SUBCASE("Bytes below 1024") {
    const auto result = unit.formatPrefix(512.0);
    CHECK(result.find("512.") != std::string::npos);
    CHECK(result.find(" B") != std::string::npos);
  }

  SUBCASE("KiB") {
    const auto result = unit.formatPrefix(2048.0);
    CHECK(result.find("2.0") != std::string::npos);
    CHECK(result.find("KiB") != std::string::npos);
  }

  SUBCASE("MiB") {
    const auto result = unit.formatPrefix(1048576.0);
    CHECK(result.find("1.0") != std::string::npos);
    CHECK(result.find("MiB") != std::string::npos);
  }

  SUBCASE("GiB") {
    const auto result = unit.formatPrefix(1073741824.0);
    CHECK(result.find("1.0") != std::string::npos);
    CHECK(result.find("GiB") != std::string::npos);
  }
}

TEST_CASE("SIUnit formatTime" * doctest::test_suite("monitoring")) {
  const SIUnit unit("s", false);

  SUBCASE("Subsecond exact") {
    const auto result = unit.formatTime(0.001);
    CHECK(result.find("ms") != std::string::npos);
  }

  SUBCASE("Seconds only") {
    const auto result = unit.formatTime(30.0);
    CHECK(result.find("30") != std::string::npos);
    CHECK(result.find("min") == std::string::npos);
  }

  SUBCASE("Minutes and seconds") {
    const auto result = unit.formatTime(90.0);
    CHECK(result.find("1 min") != std::string::npos);
  }

  SUBCASE("Hours") {
    const auto result = unit.formatTime(3661.0);
    CHECK(result.find("1 h") != std::string::npos);
    CHECK(result.find("1 min") != std::string::npos);
  }

  SUBCASE("Days") {
    const auto result = unit.formatTime(86400.0);
    CHECK(result.find("1 d") != std::string::npos);
  }

  SUBCASE("Complex time") {
    const double t = 86400.0 + 7200.0 + 180.0 + 4.5;
    const auto result = unit.formatTime(t);
    CHECK(result.find("1 d") != std::string::npos);
    CHECK(result.find("2 h") != std::string::npos);
    CHECK(result.find("3 min") != std::string::npos);
  }

  SUBCASE("Inexact mode stops at leading unit") {
    const double t = 86400.0 + 0.5;
    const auto result = unit.formatTime(t, false);
    CHECK(result.find("1 d") != std::string::npos);
  }
}

TEST_CASE("SIUnit formatScientific" * doctest::test_suite("monitoring")) {
  const SIUnit unit("J", false);

  SUBCASE("Basic") {
    const auto result = unit.formatScientific(1.5e10);
    CHECK(result.find("1.5") != std::string::npos);
    CHECK(result.find("e+") != std::string::npos);
    CHECK(result.find(" J") != std::string::npos);
  }

  SUBCASE("With error") {
    const auto result = unit.formatScientific(1.5e10, 0.1e10);
    CHECK(result.find('(') != std::string::npos);
    CHECK(result.find(')') != std::string::npos);
  }

  SUBCASE("Negative") {
    const auto result = unit.formatScientific(-3.14);
    CHECK(result.find("-3.") != std::string::npos);
    CHECK(result.find(" J") != std::string::npos);
  }
}

} // namespace seissol::unit_test

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "TestHelper.h"

#include <array>
#include <cstddef>
#include <memory>
#include <tuple>
#include <utility>

namespace seissol::unit_test {

using namespace seissol::io::datatype;

// ---------------------------------------------------------------------------
// Primitive inference
// ---------------------------------------------------------------------------

TEST_CASE("inferDatatype primitives" * doctest::test_suite("io")) {

  SUBCASE("float -> F32Datatype") {
    auto dt = inferDatatype<float>();
    CHECK(dt->size() == sizeof(float));
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "f32");
  }

  SUBCASE("double -> F64Datatype") {
    auto dt = inferDatatype<double>();
    CHECK(dt->size() == sizeof(double));
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "f64");
  }

  SUBCASE("long double -> F80Datatype") {
    auto dt = inferDatatype<long double>();
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "f80");
  }

  SUBCASE("int -> IntegerDatatype, signed") {
    auto dt = inferDatatype<int>();
    CHECK(dt->size() == sizeof(int));
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "int");
    CHECK(node["sign"].as<bool>() == true);
  }

  SUBCASE("unsigned int -> IntegerDatatype, unsigned") {
    auto dt = inferDatatype<unsigned int>();
    CHECK(dt->size() == sizeof(unsigned int));
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "int");
    CHECK(node["sign"].as<bool>() == false);
  }

  SUBCASE("char -> IntegerDatatype") {
    auto dt = inferDatatype<char>();
    CHECK(dt->size() == sizeof(char));
  }

  SUBCASE("size_t -> IntegerDatatype") {
    auto dt = inferDatatype<std::size_t>();
    CHECK(dt->size() == sizeof(std::size_t));
    const auto node = dt->serialize();
    CHECK(node["sign"].as<bool>() == false);
  }
}

// ---------------------------------------------------------------------------
// C-style array inference
// ---------------------------------------------------------------------------

TEST_CASE("inferDatatype C-arrays" * doctest::test_suite("io")) {

  SUBCASE("double[3] -> ArrayDatatype") {
    auto dt = inferDatatype<double[3]>();
    CHECK(dt->size() == 3 * sizeof(double));
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "array");
    CHECK(node["shape"].as<std::vector<std::size_t>>() == std::vector<std::size_t>{3});
  }

  SUBCASE("double[3][4] -> ArrayDatatype multidimensional") {
    auto dt = inferDatatype<double[3][4]>();
    CHECK(dt->size() == 3 * 4 * sizeof(double));
    const auto node = dt->serialize();
    CHECK(node["shape"].as<std::vector<std::size_t>>() == std::vector<std::size_t>{3, 4});
  }

  SUBCASE("int[2][3][4] -> ArrayDatatype 3D") {
    auto dt = inferDatatype<int[2][3][4]>();
    CHECK(dt->size() == 2 * 3 * 4 * sizeof(int));
    const auto node = dt->serialize();
    CHECK(node["shape"].as<std::vector<std::size_t>>() == std::vector<std::size_t>{2, 3, 4});
  }
}

// ---------------------------------------------------------------------------
// std::array inference
// ---------------------------------------------------------------------------

TEST_CASE("inferDatatype std::array" * doctest::test_suite("io")) {

  SUBCASE("std::array<float, 5>") {
    auto dt = inferDatatype<std::array<float, 5>>();
    CHECK(dt->size() == 5 * sizeof(float));
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "array");
    CHECK(node["shape"].as<std::vector<std::size_t>>() == std::vector<std::size_t>{5});
  }

  SUBCASE("std::array<double, 1>") {
    auto dt = inferDatatype<std::array<double, 1>>();
    CHECK(dt->size() == sizeof(double));
  }
}

// ---------------------------------------------------------------------------
// std::pair inference
// ---------------------------------------------------------------------------

TEST_CASE("inferDatatype std::pair" * doctest::test_suite("io")) {
  auto dt = inferDatatype<std::pair<int, double>>();
  const auto node = dt->serialize();
  CHECK(node["type"].as<std::string>() == "struct");

  // Should have two members: "first" and "second"
  const auto& members = node["members"];
  REQUIRE(members.size() == 2);
  CHECK(members[0]["name"].as<std::string>() == "first");
  CHECK(members[1]["name"].as<std::string>() == "second");
}

// ---------------------------------------------------------------------------
// std::tuple inference
// ---------------------------------------------------------------------------

TEST_CASE("inferDatatype std::tuple" * doctest::test_suite("io")) {

  SUBCASE("Two-element tuple") {
    auto dt = inferDatatype<std::tuple<int, float>>();
    const auto node = dt->serialize();
    CHECK(node["type"].as<std::string>() == "struct");
    CHECK(node["members"].size() == 2);
    CHECK(node["members"][0]["name"].as<std::string>() == "_0");
    CHECK(node["members"][1]["name"].as<std::string>() == "_1");
  }

  SUBCASE("Three-element tuple") {
    auto dt = inferDatatype<std::tuple<int, float, double>>();
    const auto node = dt->serialize();
    CHECK(node["members"].size() == 3);
    CHECK(node["members"][2]["name"].as<std::string>() == "_2");
  }
}

// ---------------------------------------------------------------------------
// Opaque fallback for non-matching types
// ---------------------------------------------------------------------------

TEST_CASE("inferDatatype opaque fallback" * doctest::test_suite("io")) {
  struct CustomPOD {
    char a;
    char b;
    char c;
  };

  auto dt = inferDatatype<CustomPOD>();
  CHECK(dt->size() == sizeof(CustomPOD));
  const auto node = dt->serialize();
  CHECK(node["type"].as<std::string>() == "opaque");
}

} // namespace seissol::unit_test

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Datatype.h"
#include "TestHelper.h"

#include <array>
#include <cstddef>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace seissol::unit_test {

using namespace seissol::io::datatype;

// ---------------------------------------------------------------------------
// Primitive types: F32, F64, F80
// ---------------------------------------------------------------------------

TEST_CASE("F32Datatype" * doctest::test_suite("io")) {
  F32Datatype dt;

  SUBCASE("Size") { CHECK(dt.size() == 4); }

  SUBCASE("Serialization") {
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "f32");
  }

  SUBCASE("toString roundtrip") {
    float val = 3.14f;
    const auto str = dt.toStringRaw(&val);
    CHECK(!str.empty());

    const auto back = dt.fromStringRaw(str);
    REQUIRE(back.has_value());
    CHECK(back->size() == sizeof(float));

    float result = 0;
    std::memcpy(&result, back->data(), sizeof(float));
    CHECK(result == doctest::Approx(val).epsilon(1e-6));
  }
}

TEST_CASE("F64Datatype" * doctest::test_suite("io")) {
  F64Datatype dt;

  SUBCASE("Size") { CHECK(dt.size() == 8); }

  SUBCASE("Serialization") {
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "f64");
  }

  SUBCASE("toString roundtrip") {
    double val = 2.718281828459045;
    const auto str = dt.toStringRaw(&val);
    const auto back = dt.fromStringRaw(str);
    REQUIRE(back.has_value());

    double result = 0;
    std::memcpy(&result, back->data(), sizeof(double));
    CHECK(result == doctest::Approx(val).epsilon(1e-14));
  }
}

TEST_CASE("F80Datatype" * doctest::test_suite("io")) {
  F80Datatype dt;

  SUBCASE("Size") { CHECK(dt.size() == 10); }

  SUBCASE("Serialization") {
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "f80");
  }
}

// ---------------------------------------------------------------------------
// IntegerDatatype
// ---------------------------------------------------------------------------

TEST_CASE("IntegerDatatype" * doctest::test_suite("io")) {

  SUBCASE("Signed 4 byte") {
    IntegerDatatype dt(4, true);
    CHECK(dt.size() == 4);
    CHECK(dt.sign() == true);
  }

  SUBCASE("Unsigned 8 byte") {
    IntegerDatatype dt(8, false);
    CHECK(dt.size() == 8);
    CHECK(dt.sign() == false);
  }

  SUBCASE("Serialization roundtrip") {
    IntegerDatatype dt(4, true);
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "int");
    CHECK(node["size"].as<std::size_t>() == 4);
    CHECK(node["sign"].as<bool>() == true);

    IntegerDatatype dt2(node);
    CHECK(dt2.size() == 4);
    CHECK(dt2.sign() == true);
  }

  SUBCASE("toString roundtrip") {
    IntegerDatatype dt(sizeof(long long), true);
    long long val = 42;
    const auto str = dt.toStringRaw(&val);
    CHECK(str == "42");

    const auto back = dt.fromStringRaw("42");
    REQUIRE(back.has_value());

    long long result = 0;
    std::memcpy(&result, back->data(), sizeof(long long));
    CHECK(result == 42);
  }
}

// ---------------------------------------------------------------------------
// StringDatatype
// ---------------------------------------------------------------------------

TEST_CASE("StringDatatype" * doctest::test_suite("io")) {
  StringDatatype dt(5);

  SUBCASE("Size") { CHECK(dt.size() == 5); }

  SUBCASE("Serialization") {
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "string");
    CHECK(node["size"].as<std::size_t>() == 5);
  }

  SUBCASE("toString") {
    const char data[] = "Hello";
    const auto str = dt.toStringRaw(data);
    CHECK(str == "Hello");
  }

  SUBCASE("fromString") {
    const auto back = dt.fromStringRaw("World");
    REQUIRE(back.has_value());
    CHECK(std::string(back->begin(), back->end()) == "World");
  }
}

// ---------------------------------------------------------------------------
// OpaqueDatatype
// ---------------------------------------------------------------------------

TEST_CASE("OpaqueDatatype" * doctest::test_suite("io")) {

  SUBCASE("Size") {
    OpaqueDatatype dt(16);
    CHECK(dt.size() == 16);
  }

  SUBCASE("Serialization roundtrip") {
    OpaqueDatatype dt(7);
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "opaque");
    CHECK(node["size"].as<std::size_t>() == 7);

    OpaqueDatatype dt2(node);
    CHECK(dt2.size() == 7);
  }

  SUBCASE("Base64 roundtrip 3 bytes") {
    OpaqueDatatype dt(3);
    const uint8_t data[] = {1, 2, 3};
    const auto encoded = dt.toStringRaw(data);
    CHECK(!encoded.empty());

    const auto decoded = dt.fromStringRaw(encoded);
    REQUIRE(decoded.has_value());
    CHECK(decoded->size() == 3);
    CHECK((*decoded)[0] == '\x01');
    CHECK((*decoded)[1] == '\x02');
    CHECK((*decoded)[2] == '\x03');
  }

  SUBCASE("Base64 roundtrip 6 bytes — no padding") {
    OpaqueDatatype dt(6);
    const uint8_t data[] = {0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
    const auto encoded = dt.toStringRaw(data);
    const auto decoded = dt.fromStringRaw(encoded);
    REQUIRE(decoded.has_value());
    CHECK(decoded->size() == 6);
    for (int i = 0; i < 6; ++i) {
      CHECK((*decoded)[i] == data[i]);
    }
  }
}

// ---------------------------------------------------------------------------
// ArrayDatatype
// ---------------------------------------------------------------------------

TEST_CASE("ArrayDatatype" * doctest::test_suite("io")) {
  auto base = std::make_shared<F64Datatype>();
  ArrayDatatype dt(base, {3, 4});

  SUBCASE("Size") {
    // 8 bytes * 3 * 4 = 96
    CHECK(dt.size() == 96);
  }

  SUBCASE("Dimensions") {
    const auto& dims = dt.dimensions();
    REQUIRE(dims.size() == 2);
    CHECK(dims[0] == 3);
    CHECK(dims[1] == 4);
  }

  SUBCASE("Base type") { CHECK(dt.base()->size() == 8); }

  SUBCASE("Serialization roundtrip") {
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "array");
    CHECK(node["shape"].as<std::vector<std::size_t>>() == std::vector<std::size_t>{3, 4});

    ArrayDatatype dt2(node);
    CHECK(dt2.size() == 96);
    CHECK(dt2.dimensions() == std::vector<std::size_t>{3, 4});
  }
}

// ---------------------------------------------------------------------------
// StructDatatype
// ---------------------------------------------------------------------------

TEST_CASE("StructDatatype" * doctest::test_suite("io")) {
  auto intType = std::make_shared<IntegerDatatype>(4, true);
  auto floatType = std::make_shared<F32Datatype>();

  std::vector<StructDatatype::MemberInfo> members = {
      {"x", 0, intType},
      {"y", 4, floatType},
  };

  SUBCASE("Auto-computed minSize") {
    StructDatatype dt(members);
    // minSize = max(0+4, 4+4) = 8
    CHECK(dt.size() == 8);
  }

  SUBCASE("Explicit size >= minSize") {
    StructDatatype dt(members, 16);
    CHECK(dt.size() == 16);
  }

  SUBCASE("Members accessor") {
    StructDatatype dt(members);
    const auto& m = dt.members();
    REQUIRE(m.size() == 2);
    CHECK(m[0].name == "x");
    CHECK(m[0].offset == 0);
    CHECK(m[1].name == "y");
    CHECK(m[1].offset == 4);
  }

  SUBCASE("Serialization roundtrip") {
    StructDatatype dt(members, 12);
    const auto node = dt.serialize();
    CHECK(node["type"].as<std::string>() == "struct");
    CHECK(node["size"].as<std::size_t>() == 12);

    StructDatatype dt2(node);
    CHECK(dt2.size() == 12);
    CHECK(dt2.members().size() == 2);
    CHECK(dt2.members()[0].name == "x");
    CHECK(dt2.members()[1].name == "y");
  }
}

// ---------------------------------------------------------------------------
// Datatype::deserialize dispatch
// ---------------------------------------------------------------------------

TEST_CASE("Datatype deserialize dispatch" * doctest::test_suite("io")) {

  SUBCASE("f32") {
    YAML::Node node;
    node["type"] = "f32";
    auto dt = Datatype::deserialize(node);
    CHECK(dt->size() == 4);
  }

  SUBCASE("f64") {
    YAML::Node node;
    node["type"] = "f64";
    auto dt = Datatype::deserialize(node);
    CHECK(dt->size() == 8);
  }

  SUBCASE("int") {
    YAML::Node node;
    node["type"] = "int";
    node["size"] = 4;
    node["sign"] = true;
    auto dt = Datatype::deserialize(node);
    CHECK(dt->size() == 4);
  }

  SUBCASE("opaque") {
    YAML::Node node;
    node["type"] = "opaque";
    node["size"] = 16;
    auto dt = Datatype::deserialize(node);
    CHECK(dt->size() == 16);
  }

  SUBCASE("Unknown type falls back to OpaqueDatatype(0)") {
    YAML::Node node;
    node["type"] = "unknown_type";
    auto dt = Datatype::deserialize(node);
    CHECK(dt->size() == 0);
  }
}

} // namespace seissol::unit_test

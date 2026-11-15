// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Filesystem.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <fstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {
using namespace seissol::initializer::parameters;

// ---------------------------------------------------------------------------
// sanitize()
// ---------------------------------------------------------------------------

TEST_CASE("sanitize" * doctest::test_suite("initializer")) {
  SUBCASE("Lowercase") {
    std::string s = "HELLO";
    sanitize(s);
    CHECK(s == "hello");
  }

  SUBCASE("Trim whitespace") {
    std::string s = "  world  ";
    sanitize(s);
    CHECK(s == "world");
  }

  SUBCASE("Trim and lowercase combined") {
    std::string s = "  FoO BaR  ";
    sanitize(s);
    CHECK(s == "foo bar");
  }

  SUBCASE("Already clean") {
    std::string s = "clean";
    sanitize(s);
    CHECK(s == "clean");
  }

  SUBCASE("Empty string") {
    std::string s;
    sanitize(s);
    CHECK(s.empty());
  }
}

// ---------------------------------------------------------------------------
// ParameterReader: basic read operations
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader read<T>" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    intfield: 42
    doublefield: 3.14
    strfield: hello
    boolfield: 1
    boolfalse: 0
  )");
  ParameterReader reader(node, "", false);

  SUBCASE("Read int") {
    auto val = reader.read<int>("intfield");
    REQUIRE(val.has_value());
    CHECK(val.value() == 42);
  }

  SUBCASE("Read double") {
    auto val = reader.read<double>("doublefield");
    REQUIRE(val.has_value());
    CHECK(val.value() == doctest::Approx(3.14));
  }

  SUBCASE("Read string") {
    auto val = reader.read<std::string>("strfield");
    REQUIRE(val.has_value());
    CHECK(val.value() == "hello");
  }

  SUBCASE("Read bool true") {
    auto val = reader.read<bool>("boolfield");
    REQUIRE(val.has_value());
    CHECK(val.value() == true);
  }

  SUBCASE("Read bool false") {
    auto val = reader.read<bool>("boolfalse");
    REQUIRE(val.has_value());
    CHECK(val.value() == false);
  }

  SUBCASE("Missing field returns empty optional") {
    const auto val = reader.read<int>("nonexistent");
    CHECK_FALSE(val.has_value());
  }
}

// ---------------------------------------------------------------------------
// ParameterReader: readWithDefault
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader readWithDefault" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load("present: 99");
  ParameterReader reader(node, "", false);

  SUBCASE("Returns value when present") { CHECK(reader.readWithDefault("present", 0) == 99); }

  SUBCASE("Returns default when missing") { CHECK(reader.readWithDefault("absent", 77) == 77); }

  SUBCASE("Default string") {
    CHECK(reader.readWithDefault<std::string>("absent", "fallback") == "fallback");
  }
}

// ---------------------------------------------------------------------------
// ParameterReader: readWithDefaultStringEnum
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader readWithDefaultStringEnum" * doctest::test_suite("initializer")) {
  enum class Color { Red, Green, Blue };
  const std::unordered_map<std::string, Color> colorMap = {
      {"red", Color::Red}, {"green", Color::Green}, {"blue", Color::Blue}};

  SUBCASE("Reads present value") {
    const YAML::Node node = YAML::Load("color: Green");
    ParameterReader reader(node, "", false);
    auto result = reader.readWithDefaultStringEnum("color", "red", colorMap);
    CHECK(result == Color::Green);
  }

  SUBCASE("Uses default when missing") {
    const YAML::Node node = YAML::Load("other: 1");
    ParameterReader reader(node, "", false);
    auto result = reader.readWithDefaultStringEnum("color", "blue", colorMap);
    CHECK(result == Color::Blue);
  }

  SUBCASE("Case-insensitive via sanitize") {
    const YAML::Node node = YAML::Load("color: RED");
    ParameterReader reader(node, "", false);
    auto result = reader.readWithDefaultStringEnum("color", "blue", colorMap);
    CHECK(result == Color::Red);
  }
}

// ---------------------------------------------------------------------------
// ParameterReader: empty reader
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader empty mode" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load("field: 42");
  ParameterReader reader(node, "", true); // empty=true → hasField always false

  SUBCASE("Read returns empty even if data exists in node") {
    const auto val = reader.read<int>("field");
    CHECK_FALSE(val.has_value());
  }

  SUBCASE("readWithDefault always returns default") {
    CHECK(reader.readWithDefault("field", 0) == 0);
  }
}

// ---------------------------------------------------------------------------
// ParameterReader: readSubNode
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader readSubNode" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    equations:
      materialfilename: rock.yaml
      plasticity: 0
  )");
  ParameterReader reader(node, "", false);

  auto* sub = reader.readSubNode("equations");
  REQUIRE(sub != nullptr);

  SUBCASE("Read from subnode") {
    auto mat = sub->read<std::string>("materialfilename");
    REQUIRE(mat.has_value());
    CHECK(mat.value() == "rock.yaml");
  }

  SUBCASE("Read bool from subnode") {
    auto val = sub->read<bool>("plasticity");
    REQUIRE(val.has_value());
    CHECK(val.value() == false);
  }

  SUBCASE("Missing subnode returns empty reader") {
    auto* missing = reader.readSubNode("nonexistentsection");
    REQUIRE(missing != nullptr);
    const auto val = missing->read<int>("field");
    CHECK_FALSE(val.has_value());
  }
}

// ---------------------------------------------------------------------------
// ParameterReader: hasField
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader hasField" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load("alpha: 1\nbeta: 2");
  ParameterReader reader(node, "", false);

  CHECK(reader.hasField("alpha"));
  CHECK(reader.hasField("beta"));
  CHECK_FALSE(reader.hasField("gamma"));
}

// ---------------------------------------------------------------------------
// ParameterReader: readPath
// ---------------------------------------------------------------------------

TEST_CASE("ParameterReader readPath" * doctest::test_suite("initializer")) {
  const std::string directory = "readpath-test";
  std::filesystem::create_directories(directory);
  {
    std::ofstream file(directory + "/material.yaml");
    file << "!ConstantMap\n";
  }
  const YAML::Node node = YAML::Load(R"(
    material: material.yaml
    missing: absent.yaml
    empty: ''
  )");
  // a parameter file in a directory of its own, as a job script names it
  ParameterReader reader(node, directory + "/parameters.par", false);

  SUBCASE("A file next to the parameter file") {
    const auto path = reader.readPath("material");
    REQUIRE(path.has_value());
    CHECK(path.value() == std::filesystem::canonical(directory + "/material.yaml").string());
  }

  SUBCASE("A file that does not exist comes back as given") {
    CHECK(reader.readPath("missing") == std::optional<std::string>("absent.yaml"));
  }

  SUBCASE("An empty name stays empty") {
    CHECK(reader.readPath("empty") == std::optional<std::string>(""));
  }

  SUBCASE("A name that is not given") { CHECK_FALSE(reader.readPath("other").has_value()); }

  std::filesystem::remove_all(directory);
}

TEST_CASE("ParameterReader paths and unknown fields" * doctest::test_suite("initializer")) {
  const YAML::Node node = YAML::Load(R"(
    equations:
      materialfilename: rock.yaml
      typo: 1
    output:
      rfilename: receivers.dat
    unread:
      field: 1
  )");
  ParameterReader reader(node, "", false);
  reader.readSubNode("output")->readPath("rfilename");
  reader.readSubNode("equations")->readPath("materialfilename");

  const auto paths = reader.paths();
  REQUIRE(paths.size() == 2);
  // by section, in alphabetical order
  CHECK(paths[0].field == std::vector<std::string>{"equations", "materialfilename"});
  CHECK(paths[0].path == "rock.yaml");
  CHECK(paths[1].field == std::vector<std::string>{"output", "rfilename"});
  CHECK(paths[1].path == "receivers.dat");

  const auto unknown = reader.unknownFields();
  REQUIRE(unknown.size() == 2);
  CHECK(unknown[0] == "unread");
  CHECK(unknown[1] == "equations.typo");
}

} // namespace seissol::unit_test

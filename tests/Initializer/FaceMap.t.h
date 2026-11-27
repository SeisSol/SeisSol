// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/BasicTypedefs.h"
#include "Initializer/FaceMap.h"

#include <cstdio>
#include <fstream>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::unit_test {

TEST_CASE("Default face map" * doctest::test_suite("initializer")) {
  const auto map = defaultFaceMap();

  CHECK(map.at(0) == FaceType::Regular);
  CHECK(map.at(1) == FaceType::FreeSurface);
  CHECK(map.at(2) == FaceType::FreeSurfaceGravity);
  CHECK(map.at(3) == FaceType::DynamicRupture);
  CHECK(map.at(4) == FaceType::Dirichlet);
  CHECK(map.at(5) == FaceType::Outflow);
  CHECK(map.at(6) == FaceType::Regular);
  CHECK(map.at(7) == FaceType::Analytical);
  for (std::size_t i = 8; i <= 64; ++i) {
    CHECK_FALSE(map.at(i).has_value());
  }
  CHECK(map.at(65) == FaceType::DynamicRupture);
  CHECK(map.at(66) == FaceType::DynamicRupture);
  CHECK(map.at(100) == FaceType::DynamicRupture);
  CHECK(map.at(1000) == FaceType::DynamicRupture);
  CHECK(map.at(10000) == FaceType::DynamicRupture);
}

TEST_CASE("Default face map as input" * doctest::test_suite("initializer")) {
  const auto node = YAML::Load(R"(
        regular:
            - 0
            - 6
        freeSurface:
            - 1
        freeSurfaceGravity:
            - 2
        dynamicRupture:
            - 3
            - 65- # note the open range
        dirichlet:
            - 4
        outflow:
            - 5
        analytical:
            - 7
    )");
  const auto map = parseFaceMap(node);

  CHECK(map.at(0) == FaceType::Regular);
  CHECK(map.at(1) == FaceType::FreeSurface);
  CHECK(map.at(2) == FaceType::FreeSurfaceGravity);
  CHECK(map.at(3) == FaceType::DynamicRupture);
  CHECK(map.at(4) == FaceType::Dirichlet);
  CHECK(map.at(5) == FaceType::Outflow);
  CHECK(map.at(6) == FaceType::Regular);
  CHECK(map.at(7) == FaceType::Analytical);
  for (std::size_t i = 8; i <= 64; ++i) {
    CHECK_FALSE(map.at(i).has_value());
  }
  CHECK(map.at(65) == FaceType::DynamicRupture);
  CHECK(map.at(66) == FaceType::DynamicRupture);
  CHECK(map.at(100) == FaceType::DynamicRupture);
  CHECK(map.at(1000) == FaceType::DynamicRupture);
  CHECK(map.at(10000) == FaceType::DynamicRupture);
}

TEST_CASE("Custom face map input" * doctest::test_suite("initializer")) {
  const auto node = YAML::Load(R"(
        regular: 1, 3
        dynamicRupture:
            - 4
            - 6-10
            - 102-
        outflow: [2,5,11]
        dirichlet: 23
        freeSurface: 40-41, 44
        freeSurfaceGravity: [50-51, 54]
        analytical: []
    )");
  const auto map = parseFaceMap(node);

  CHECK_FALSE(map.at(0).has_value());
  CHECK_FALSE(map.at(12).has_value());
  CHECK_FALSE(map.at(13).has_value());

  CHECK(map.at(1) == FaceType::Regular);
  CHECK(map.at(3) == FaceType::Regular);

  CHECK(map.at(4) == FaceType::DynamicRupture);
  CHECK(map.at(6) == FaceType::DynamicRupture);
  CHECK(map.at(7) == FaceType::DynamicRupture);
  CHECK(map.at(8) == FaceType::DynamicRupture);
  CHECK(map.at(9) == FaceType::DynamicRupture);
  CHECK(map.at(10) == FaceType::DynamicRupture);
  CHECK(map.at(102) == FaceType::DynamicRupture);
  CHECK(map.at(103) == FaceType::DynamicRupture);

  CHECK(map.at(2) == FaceType::Outflow);
  CHECK(map.at(5) == FaceType::Outflow);
  CHECK(map.at(11) == FaceType::Outflow);

  CHECK(map.at(23) == FaceType::Dirichlet);

  CHECK(map.at(40) == FaceType::FreeSurface);
  CHECK(map.at(41) == FaceType::FreeSurface);
  CHECK(map.at(44) == FaceType::FreeSurface);

  CHECK(map.at(50) == FaceType::FreeSurfaceGravity);
  CHECK(map.at(51) == FaceType::FreeSurfaceGravity);
  CHECK(map.at(54) == FaceType::FreeSurfaceGravity);
}

TEST_CASE("Face map from a file" * doctest::test_suite("initializer")) {
  const std::string path = "facemap-test.yaml";
  {
    std::ofstream file(path);
    file << "regular: 0\nfreeSurface: 9\ndynamicRupture: 3, 65-\n";
  }
  const auto map = readFaceMap(path);
  CHECK(std::remove(path.c_str()) == 0);

  CHECK(map.at(0) == FaceType::Regular);
  CHECK(map.at(9) == FaceType::FreeSurface);
  CHECK(map.at(3) == FaceType::DynamicRupture);
  CHECK(map.at(70) == FaceType::DynamicRupture);
  CHECK_FALSE(map.at(1).has_value());
}

TEST_CASE("Face map entries with a flux" * doctest::test_suite("initializer")) {
  using initializer::parameters::FluxSetting;
  using initializer::parameters::NumericalFlux;
  const auto node = YAML::Load(R"(
        regular:
            - 0, 6
            - tags: 100-110
              numflux: godunov-shearcentered
            - groups: 2, 5-
              numfluxnearfault: mixed
              fluxgpartnearfault: 0.5
              fluxdampingnearfault: 1 1 1 0.5
            - tags: [120, 130-139]
            - tags: 140
              groups: 3
              fluxdamping: 1 0 1 1
        dynamicRupture:
            - 3
            - tags: 65-99
    )");

  // the tags of the entries that are maps belong to their face type
  const auto map = parseFaceMap(node);
  CHECK(map.at(0) == FaceType::Regular);
  CHECK(map.at(100) == FaceType::Regular);
  CHECK(map.at(110) == FaceType::Regular);
  CHECK(map.at(120) == FaceType::Regular);
  CHECK(map.at(135) == FaceType::Regular);
  CHECK(map.at(140) == FaceType::Regular);
  CHECK(map.at(70) == FaceType::DynamicRupture);
  CHECK_FALSE(map.at(111).has_value());

  // the entries that give a flux are its rules, in their order; tags alone are none
  const auto rules = parseFaceFluxRules(node);
  REQUIRE(rules.size() == 3);

  CHECK(rules[0].tags.has_value());
  CHECK_FALSE(rules[0].groups.has_value());
  CHECK(rules[0].flux ==
        FluxSetting{
            .flux = NumericalFlux::Mixed, .gpart = 1, .rpart = 0.1, .damping = {1, 1, 1, 0, 1, 1}});
  CHECK_FALSE(rules[0].fluxNearFault.has_value());
  CHECK(rules[0].selects(100, 0, 0));
  CHECK(rules[0].selects(110, 7, 8));
  CHECK_FALSE(rules[0].selects(111, 0, 0));

  CHECK_FALSE(rules[1].tags.has_value());
  CHECK_FALSE(rules[1].flux.has_value());
  CHECK(rules[1].fluxNearFault == FluxSetting{.flux = NumericalFlux::Mixed,
                                              .gpart = 0.5,
                                              .rpart = 0,
                                              .damping = {1, 1, 1, 0.5, 1, 1}});
  // a cell of one of the groups on either side
  CHECK(rules[1].selects(0, 2, 1));
  CHECK(rules[1].selects(0, 1, 7));
  CHECK_FALSE(rules[1].selects(0, 1, 3));

  CHECK(rules[2].flux ==
        FluxSetting{.flux = NumericalFlux::Godunov, .damping = {1, 0, 1, 1, 1, 1}});
  CHECK(rules[2].selects(140, 3, 1));
  CHECK_FALSE(rules[2].selects(140, 1, 1));
  CHECK_FALSE(rules[2].selects(0, 3, 3));
}

} // namespace seissol::unit_test

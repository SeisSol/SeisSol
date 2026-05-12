// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Initializer/BasicTypedefs.h"
#include "Initializer/FaceMap.h"

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

} // namespace seissol::unit_test

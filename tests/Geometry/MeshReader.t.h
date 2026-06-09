// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Constants.h"
#include "Geometry/MeshReader.h"
#include "TestHelper.h"

#include <Eigen/Dense>
#include <array>

namespace seissol::unit_test {
class MockReader2 : public seissol::geometry::MeshReader {
  public:
  explicit MockReader2(FaceType boundaryType) {
    vertices_.resize(4);

    elements_.resize(1);
    elements_.at(0).vertices[0] = 0;
    elements_.at(0).vertices[1] = 1;
    elements_.at(0).vertices[2] = 2;
    elements_.at(0).vertices[3] = 3;
    elements_.at(0).boundaries[0] = boundaryType;
    elements_.at(0).boundaries[1] = boundaryType;
    elements_.at(0).boundaries[2] = boundaryType;
    elements_.at(0).boundaries[3] = boundaryType;
  }
};

TEST_CASE("MeshReader" * doctest::test_suite("geometry")) {
  SUBCASE("No DR") {
    MockReader2 rdr(FaceType::DynamicRupture);

    rdr.disableDR();

    REQUIRE(rdr.getElements()[0].boundaries[0] == FaceType::Regular);
    REQUIRE(rdr.getElements()[0].boundaries[1] == FaceType::Regular);
    REQUIRE(rdr.getElements()[0].boundaries[2] == FaceType::Regular);
    REQUIRE(rdr.getElements()[0].boundaries[3] == FaceType::Regular);
  }
}

} // namespace seissol::unit_test

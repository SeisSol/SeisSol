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
  explicit MockReader2(int boundaryType) : seissol::geometry::MeshReader(0) {
    m_vertices.resize(4);

    m_elements.resize(1);
    m_elements.at(0).vertices[0] = 0;
    m_elements.at(0).vertices[1] = 1;
    m_elements.at(0).vertices[2] = 2;
    m_elements.at(0).vertices[3] = 3;
    m_elements.at(0).boundaries[0] = boundaryType;
    m_elements.at(0).boundaries[1] = boundaryType;
    m_elements.at(0).boundaries[2] = boundaryType;
    m_elements.at(0).boundaries[3] = boundaryType;
  }
};

TEST_CASE("MeshReader") {
  SUBCASE("No DR") {
    MockReader2 rdr(3);

    rdr.disableDR();

    REQUIRE(rdr.getElements()[0].boundaries[0] == 0);
    REQUIRE(rdr.getElements()[0].boundaries[1] == 0);
    REQUIRE(rdr.getElements()[0].boundaries[2] == 0);
    REQUIRE(rdr.getElements()[0].boundaries[3] == 0);
  }
}

} // namespace seissol::unit_test

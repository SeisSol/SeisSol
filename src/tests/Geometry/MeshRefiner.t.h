// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "tests/TestHelper.h"
#include <Eigen/Dense>
#include <array>

#include "Geometry/Refinement/MeshRefiner.h"
#include "Geometry/Refinement/RefinerUtils.h"
#include "MockReader.h"

namespace seissol::unit_test {

inline void assertPoint(const double* a, const Eigen::Vector3d& b, double epsilon) {
  for (int i = 0; i < 3; i++) {
    REQUIRE(a[i] == AbsApprox(b[i]).epsilon(epsilon));
  }
}

inline void assertCell(const unsigned int* a, const Eigen::Vector4i& b) {
  for (int i = 0; i < 4; i++) {
    REQUIRE(a[i] == b[i]);
  }
}

TEST_CASE("Mesh refiner") {
  constexpr double Epsilon = std::numeric_limits<double>::epsilon();
  std::srand(0);
  const auto vertices =
      std::array<Eigen::Vector3d, 4>{{Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX)}};
  const seissol::MockReader mockReader(vertices);

  SUBCASE("Divide by 4") {
    const std::array<Eigen::Vector3d, 5> expectedVerticesDivideBy4 = {
        vertices[0],
        vertices[1],
        vertices[2],
        vertices[3],
        0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3])};

    const std::array<Eigen::Vector4i, 4> expectedCellsDivideBy4{Eigen::Vector4i(0, 1, 2, 4),
                                                                Eigen::Vector4i(0, 1, 3, 4),
                                                                Eigen::Vector4i(0, 2, 3, 4),
                                                                Eigen::Vector4i(1, 2, 3, 4)};

    const seissol::refinement::DivideTetrahedronBy4<double> refineBy4;
    const seissol::refinement::MeshRefiner<double> meshRefiner(mockReader, refineBy4);
    REQUIRE(meshRefiner.getNumCells() == 4);
    REQUIRE(meshRefiner.getNumVertices() == 5);
    for (unsigned i = 0; i < meshRefiner.getNumVertices(); i++) {
      assertPoint(&meshRefiner.getVertexData()[3 * i], expectedVerticesDivideBy4[i], Epsilon);
    }
    for (unsigned i = 0; i < meshRefiner.getNumCells(); i++) {
      assertCell(&meshRefiner.getCellData()[4 * i], expectedCellsDivideBy4[i]);
    }
  }
  SUBCASE("Divide by 8") {
    const std::array<Eigen::Vector3d, 10> expectedVerticesDivideBy8 = {
        vertices[0],
        vertices[1],
        vertices[2],
        vertices[3],
        0.5 * (vertices[0] + vertices[1]),
        0.5 * (vertices[0] + vertices[2]),
        0.5 * (vertices[0] + vertices[3]),
        0.5 * (vertices[1] + vertices[2]),
        0.5 * (vertices[1] + vertices[3]),
        0.5 * (vertices[2] + vertices[3])};

    const std::array<Eigen::Vector4i, 8> expectedCellsDivideBy8{Eigen::Vector4i(0, 4, 5, 6),
                                                                Eigen::Vector4i(1, 4, 7, 8),
                                                                Eigen::Vector4i(2, 5, 7, 9),
                                                                Eigen::Vector4i(3, 6, 8, 9),
                                                                Eigen::Vector4i(4, 5, 6, 8),
                                                                Eigen::Vector4i(4, 5, 7, 8),
                                                                Eigen::Vector4i(5, 6, 8, 9),
                                                                Eigen::Vector4i(5, 7, 8, 9)};

    const seissol::refinement::DivideTetrahedronBy8<double> refineBy8;
    const seissol::refinement::MeshRefiner<double> meshRefiner(mockReader, refineBy8);
    REQUIRE(meshRefiner.getNumCells() == 8);
    REQUIRE(meshRefiner.getNumVertices() == 10);
    for (unsigned i = 0; i < meshRefiner.getNumVertices(); i++) {
      assertPoint(&meshRefiner.getVertexData()[3 * i], expectedVerticesDivideBy8[i], Epsilon);
    }
    for (unsigned i = 0; i < meshRefiner.getNumCells(); i++) {
      assertCell(&meshRefiner.getCellData()[4 * i], expectedCellsDivideBy8[i]);
    }
  }

  SUBCASE("Divide by 32") {
    const std::array<Eigen::Vector3d, 18> expectedVerticesDivideBy32 = {
        vertices[0],
        vertices[1],
        vertices[2],
        vertices[3],
        0.5 * (vertices[0] + vertices[1]),
        0.5 * (vertices[0] + vertices[2]),
        0.5 * (vertices[0] + vertices[3]),
        0.5 * (vertices[1] + vertices[2]),
        0.5 * (vertices[1] + vertices[3]),
        0.5 * (vertices[2] + vertices[3]),

        0.625 * vertices[0] + 0.125 * vertices[1] + 0.125 * vertices[2] + 0.125 * vertices[3],
        0.125 * vertices[0] + 0.625 * vertices[1] + 0.125 * vertices[2] + 0.125 * vertices[3],
        0.125 * vertices[0] + 0.125 * vertices[1] + 0.625 * vertices[2] + 0.125 * vertices[3],
        0.125 * vertices[0] + 0.125 * vertices[1] + 0.125 * vertices[2] + 0.625 * vertices[3],

        0.375 * vertices[0] + 0.25 * vertices[1] + 0.125 * vertices[2] + 0.25 * vertices[3],
        0.25 * vertices[0] + 0.375 * vertices[1] + 0.25 * vertices[2] + 0.125 * vertices[3],
        0.25 * vertices[0] + 0.125 * vertices[1] + 0.25 * vertices[2] + 0.375 * vertices[3],
        0.125 * vertices[0] + 0.25 * vertices[1] + 0.375 * vertices[2] + 0.25 * vertices[3]};

    const std::array<Eigen::Vector4i, 32> expectedCellsDivideBy32{
        Eigen::Vector4i(0, 4, 5, 10), Eigen::Vector4i(0, 4, 6, 10), Eigen::Vector4i(0, 5, 6, 10),
        Eigen::Vector4i(4, 5, 6, 10), Eigen::Vector4i(1, 4, 7, 11), Eigen::Vector4i(1, 4, 8, 11),
        Eigen::Vector4i(1, 7, 8, 11), Eigen::Vector4i(4, 7, 8, 11), Eigen::Vector4i(2, 5, 7, 12),
        Eigen::Vector4i(2, 5, 9, 12), Eigen::Vector4i(2, 7, 9, 12), Eigen::Vector4i(5, 7, 9, 12),
        Eigen::Vector4i(3, 6, 8, 13), Eigen::Vector4i(3, 6, 9, 13), Eigen::Vector4i(3, 8, 9, 13),
        Eigen::Vector4i(6, 8, 9, 13), Eigen::Vector4i(4, 5, 6, 14), Eigen::Vector4i(4, 5, 8, 14),
        Eigen::Vector4i(4, 6, 8, 14), Eigen::Vector4i(5, 6, 8, 14), Eigen::Vector4i(4, 5, 7, 15),
        Eigen::Vector4i(4, 5, 8, 15), Eigen::Vector4i(4, 7, 8, 15), Eigen::Vector4i(5, 7, 8, 15),
        Eigen::Vector4i(5, 6, 8, 16), Eigen::Vector4i(5, 6, 9, 16), Eigen::Vector4i(5, 8, 9, 16),
        Eigen::Vector4i(6, 8, 9, 16), Eigen::Vector4i(5, 7, 8, 17), Eigen::Vector4i(5, 7, 9, 17),
        Eigen::Vector4i(5, 8, 9, 17), Eigen::Vector4i(7, 8, 9, 17)};

    const seissol::refinement::DivideTetrahedronBy32<double> refineBy32;
    const seissol::refinement::MeshRefiner<double> meshRefiner(mockReader, refineBy32);
    REQUIRE(meshRefiner.getNumCells() == 32);
    REQUIRE(meshRefiner.getNumVertices() == 18);
    for (unsigned i = 0; i < meshRefiner.getNumVertices(); i++) {
      assertPoint(&meshRefiner.getVertexData()[3 * i], expectedVerticesDivideBy32[i], Epsilon);
    }
    for (unsigned i = 0; i < meshRefiner.getNumCells(); i++) {
      assertCell(&meshRefiner.getCellData()[4 * i], expectedCellsDivideBy32[i]);
    }
  }
}

} // namespace seissol::unit_test

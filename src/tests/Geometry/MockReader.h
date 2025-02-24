// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_TESTS_GEOMETRY_MOCKREADER_H_
#define SEISSOL_SRC_TESTS_GEOMETRY_MOCKREADER_H_

#include <Eigen/Dense>
#include <array>

#include "Geometry/MeshReader.h"

namespace seissol {
class MockReader : public seissol::geometry::MeshReader {
  public:
  MockReader(std::array<Eigen::Vector3d, 4> vertices) : seissol::geometry::MeshReader(0) {
    m_vertices.resize(4);
    for (int i = 0; i < 4; i++) {
      std::copy(vertices[i].data(), vertices[i].data() + 3, m_vertices.at(i).coords);
      m_vertices.at(i).elements = {1};
    }

    m_elements.resize(1);
    m_elements.at(0).vertices[0] = 0;
    m_elements.at(0).vertices[1] = 1;
    m_elements.at(0).vertices[2] = 2;
    m_elements.at(0).vertices[3] = 3;
  }
};
} // namespace seissol

#endif // SEISSOL_SRC_TESTS_GEOMETRY_MOCKREADER_H_

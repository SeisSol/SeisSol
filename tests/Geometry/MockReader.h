// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_TESTS_GEOMETRY_MOCKREADER_H_
#define SEISSOL_TESTS_GEOMETRY_MOCKREADER_H_

#include "Common/Constants.h"
#include "Geometry/MeshReader.h"

#include <Eigen/Dense>
#include <array>

namespace seissol {
class MockReader : public seissol::geometry::MeshReader {
  public:
  explicit MockReader(const std::array<Eigen::Vector3d, 4>& vertices)
      : seissol::geometry::MeshReader(0) {
    m_vertices.resize(4);
    for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
      std::copy_n(vertices[i].data(), 3, m_vertices.at(i).coords.begin());
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

#endif // SEISSOL_TESTS_GEOMETRY_MOCKREADER_H_

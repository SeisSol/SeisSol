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
  explicit MockReader(const std::array<Eigen::Vector3d, 4>& vertices) {
    vertices_.resize(4);
    for (std::size_t i = 0; i < Cell::NumVertices; ++i) {
      std::copy(vertices[i].data(), vertices[i].data() + 3, vertices_.at(i).coords);
      vertices_.at(i).elements = {1};
    }

    elements_.resize(1);
    elements_.at(0).vertices[0] = 0;
    elements_.at(0).vertices[1] = 1;
    elements_.at(0).vertices[2] = 2;
    elements_.at(0).vertices[3] = 3;
  }
};
} // namespace seissol

#endif // SEISSOL_TESTS_GEOMETRY_MOCKREADER_H_

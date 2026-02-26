// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_

#include <array>
#include <vector>

namespace seissol::io::instance::geometry {

extern const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine4;
extern const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine8;
extern const std::vector<std::vector<std::array<double, 2>>> TriangleRefine4;

// subdivides a tetrahedron by the
template <std::size_t Dim>
std::vector<std::vector<std::array<double, Dim>>>
    applySubdivide(const std::vector<std::vector<std::array<double, Dim>>>& input,
                   const std::vector<std::vector<std::array<double, Dim>>>& refine) {
  std::vector<std::vector<std::array<double, Dim>>> output;
  output.reserve(input.size() * refine.size());
  for (const auto& tetrahedron : input) {
    for (const auto& refTet : refine) {
      assert(refTet.size() == Dim + 1);

      std::vector<std::array<double, Dim>> points(tetrahedron.size());

      std::array<double, Dim> offset = refTet[0];

      std::array<std::array<double, Dim>, Dim> transform{};
      for (std::size_t i = 0; i < Dim; ++i) {
        for (std::size_t j = 0; j < Dim; ++j) {
          transform[i][j] = refTet[i + 1][j] - offset[j];
        }
      }
      for (std::size_t i = 0; i < tetrahedron.size(); ++i) {
        points[i] = offset;
        for (std::size_t j = 0; j < Dim; ++j) {
          for (std::size_t k = 0; k < Dim; ++k) {
            points[i][j] += tetrahedron[i][k] * transform[k][j];
          }
        }
      }
      output.emplace_back(std::move(points));
    }
  }
  return output;
}

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_

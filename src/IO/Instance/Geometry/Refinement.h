// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_

#include "Numerical/Projection.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>

namespace seissol::io::instance::geometry {

/**
 * @brief Subdivision tables, given as the vertices of each subcell in the coordinates of the
 * reference simplex.
 *
 * The vertex order is the one of AffineMap::fromVertices, i.e. the entries of a tetrahedron are
 * the images of @f$ 0, e_1, e_2, e_3 @f$ and therefore correspond to element.vertices[0..3].
 * Every entry is oriented positively; a subdivision built from these tables consists of valid
 * (non-inverted) output cells. Note that reorienting a cell renumbers the subcells a subsequent
 * center split produces, so the tables have to stay oriented rather than being normalized after
 * the fact.
 */
extern const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine4;
extern const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine8;
extern const std::vector<std::vector<std::array<double, 2>>> TriangleRefine4;

//! @brief A single subcell, given as the affine map taking the reference simplex onto it.
template <std::size_t Dim>
using Subcell = numerical::AffineMap<Dim>;

//! @brief The (trivial) subdivision consisting of the whole reference simplex.
template <std::size_t Dim>
std::vector<Subcell<Dim>> unrefined() {
  return {Subcell<Dim>::identity()};
}

/**
 * @brief Subdivides each subcell of @p input according to @p refine .
 *
 * The subcells of the result are the images of the @p refine simplices under the respective input
 * map, i.e. the refinement is applied *inside* each existing subcell. The enumeration is
 * input-major.
 */
template <std::size_t Dim>
std::vector<Subcell<Dim>>
    subdivideMaps(const std::vector<Subcell<Dim>>& input,
                  const std::vector<std::vector<std::array<double, Dim>>>& refine) {
  std::vector<Subcell<Dim>> output;
  output.reserve(input.size() * refine.size());
  for (const auto& map : input) {
    for (const auto& simplex : refine) {
      assert(simplex.size() == Dim + 1);
      output.emplace_back(map.compose(Subcell<Dim>::fromVertices(simplex)));
    }
  }
  return output;
}

//! @brief Evaluates @p points inside each subcell.
template <std::size_t Dim>
std::vector<std::vector<std::array<double, Dim>>>
    applyMaps(const std::vector<Subcell<Dim>>& maps,
              const std::vector<std::array<double, Dim>>& points) {
  std::vector<std::vector<std::array<double, Dim>>> output;
  output.reserve(maps.size());
  for (const auto& map : maps) {
    output.emplace_back(map.apply(points));
  }
  return output;
}

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_REFINEMENT_H_

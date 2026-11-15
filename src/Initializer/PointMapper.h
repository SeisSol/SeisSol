// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_POINTMAPPER_H_
#define SEISSOL_SRC_INITIALIZER_POINTMAPPER_H_

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>

namespace seissol::initializer {

//! The tolerance the off-fault receivers are located with (see findUniqueMeshIds).
constexpr double ReceiverTolerance = 1e-3;
//! The tolerance the point sources are located with (see findUniqueMeshIds).
constexpr double SourceTolerance = 0;

/**
 * @brief The planes of the faces of a cell, by which findUniqueMeshIds decides whether a point lies
 * in the cell.
 *
 * Each face, as MeshTools numbers them, has the plane through its first vertex with its normal,
 * which is not normalized: it is twice as long as the face is large. The value of the plane
 * equation at a point is below zero on the side of the cell, zero on the plane, and above zero
 * beyond it, up to round-off.
 */
class CellPlanes {
  public:
  //! A point with a fourth coordinate of one, at which the plane equations are evaluated.
  using Point = std::array<double, Cell::Dim + 1>;

  //! How a point lies to the cell.
  struct Fit {
    //! the largest value of the plane equations at the point, i.e. that of the face whose side the
    //! point is the farthest out on; the tolerance of findUniqueMeshIds applies to it
    double value{-std::numeric_limits<double>::infinity()};
    //! the largest distance of the point from the planes, positive beyond them
    double distance{-std::numeric_limits<double>::infinity()};
  };

  //! @brief The planes of the cell with the @p vertices , in the order of the vertices of the cells
  //! of a MeshReader.
  explicit CellPlanes(const std::array<const CoordinateT*, Cell::NumVertices>& vertices);

  static Point homogeneous(const Eigen::Vector3d& point) {
    return {point(0), point(1), point(2), 1.0};
  }

  [[nodiscard]] Fit fit(const Point& point) const {
    double maxValue = -std::numeric_limits<double>::infinity();
    double maxDistance = -std::numeric_limits<double>::infinity();

#pragma omp simd reduction(max : maxValue, maxDistance)
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      double resultFace = 0;
      for (std::size_t dim = 0; dim < Cell::Dim + 1; ++dim) {
        resultFace += equations_[dim][face] * point[dim];
      }
      maxValue = std::max(maxValue, resultFace);
      maxDistance = std::max(maxDistance, resultFace / normLengths_[face]);
    }
    return {.value = maxValue, .distance = maxDistance};
  }

  /**
   * @brief What findUniqueMeshIds ranks the cell by for a point it takes the cell for: zero where
   * the point lies in the cell up to round-off, and else the value by which it misses the cell.
   */
  [[nodiscard]] double key(const Fit& fit) const {
    return fit.distance <= onFace_ ? 0.0 : fit.value;
  }

  private:
  // the coefficients of the plane equations, by dimension (the constant last) and face
  std::array<std::array<double, Cell::NumFaces>, Cell::Dim + 1> equations_{};
  std::array<double, Cell::NumFaces> normLengths_{};
  // the distance from a plane within which a point is on the plane as far as round-off can tell
  double onFace_{0};
};

/** Finds the tetrahedrons that contain the points.
 *  In "contained" we save if the point source is contained in the mesh.
 *  Returns a vector if the point source is contained in the mesh.
 *  The tolerance parameter exists to mitigate numerical errors. Making it too large should have no
 * effect on finding the right cell; but might let points outside the mesh as a whole be found to be
 * inside. However, having it too small could make some points not being assigned to cells properly.
 *  A point is taken for a cell where the CellPlanes of the cell fit it with a value of at most the
 * tolerance.
 */
std::vector<bool> findUniqueMeshIds(const Eigen::Vector3d* points,
                                    const seissol::geometry::MeshReader& mesh,
                                    std::size_t numPoints,
                                    std::size_t* meshIds,
                                    double tolerance = 0);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_POINTMAPPER_H_

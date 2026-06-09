// SPDX-FileCopyrightText: 2015 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Rettenberger
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#include "PointMapper.h"

#include "Common/Constants.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Geometry/MeshTools.h"
#include "Parallel/MPI.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <mpi.h>
#include <utility>
#include <vector>

namespace seissol::initializer {

std::vector<bool> findUniqueMeshIds(const Eigen::Vector3d* points,
                                    const seissol::geometry::MeshReader& mesh,
                                    std::size_t numPoints,
                                    std::size_t* meshIds,
                                    double tolerance) {

  const auto& vertices = mesh.getVertices();
  const auto& elements = mesh.getElements();

  auto points1 = std::vector<std::array<double, Cell::Dim + 1>>(numPoints);
  for (std::size_t point = 0; point < numPoints; ++point) {
    for (std::size_t c = 0; c < Cell::Dim; ++c) {
      points1[point][c] = points[point](c);
    }
    points1[point][Cell::Dim] = 1.0;
  }

  const auto rank = seissol::Mpi::mpi.rank();

  std::vector<std::pair<double, int>> score(
      numPoints, std::pair<double, int>(std::numeric_limits<double>::infinity(), rank));

#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < elements.size(); ++elem) {
    auto planeEquations = std::array<std::array<double, Cell::Dim + 1>, Cell::Dim + 1>();
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      VrtxCoords n{};
      VrtxCoords p{};
      MeshTools::pointOnPlane(elements[elem], face, vertices, p);
      MeshTools::normal(elements[elem], face, vertices, n);

      for (std::size_t i = 0; i < Cell::Dim; ++i) {
        planeEquations[i][face] = n[i];
      }
      planeEquations[Cell::Dim][face] = -MeshTools::dot(n, p);
    }
    for (std::size_t point = 0; point < numPoints; ++point) {
      // geometric Interpretation (up to numerical errors, hence tolerance parameter):
      // resultFace < 0: The point is inside the face (half-space).
      // resultFace = 0: The point is exactly on the face.
      // resultFace > 0: The point is outside the face.

      // we look for the face, where resultFace is the largest; i.e. the face that will the best
      // "invalidate" our membership in the cell. We also use that value as a tie breaker
      // if we find a point to be in multiple cells due to the tolerance parameter or numerical
      // inaccuracies.

      // NOLINTNEXTLINE
      double maxValue = -std::numeric_limits<double>::infinity();

#pragma omp simd reduction(max : maxValue)
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {

        double resultFace = 0;
        for (std::size_t dim = 0; dim < Cell::Dim + 1; ++dim) {
          resultFace += planeEquations[dim][face] * points1[point][dim];
        }
        maxValue = std::max(maxValue, resultFace);
      }

      if (maxValue <= tolerance) {
        // meaning: we're below tolerance to consider the cell

#pragma omp critical
        {
          // minimize the maxValue
          const auto localId = static_cast<std::size_t>(elements[elem].localId);
          if (score[point].first > maxValue) {
            score[point] = std::pair<double, int>{maxValue, rank};
            meshIds[point] = localId;
          }
        }
      }
    }
  }

  // now reduce over all ranks for the best fit (not only duplicate ranks)

  MPI_Allreduce(MPI_IN_PLACE,
                score.data(),
                score.size(),
                MPI_DOUBLE_INT,
                MPI_MINLOC,
                seissol::Mpi::mpi.comm());

  std::vector<bool> contained(numPoints);
  for (std::size_t i = 0; i < numPoints; ++i) {
    contained[i] =
        score[i].second == rank && score[i].first < std::numeric_limits<double>::infinity();
  }
  return contained;
}

} // namespace seissol::initializer

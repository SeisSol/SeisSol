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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mpi.h>
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

  // A point may lie in more than one cell: on a face or an edge the cells meet at, or within the
  // tolerance of several of them. The cell that is taken must not depend on which thread got to it
  // first, nor on how the mesh is partitioned, so the candidates are ranked by a key that is the
  // same wherever they are looked at:
  //   - a cell that holds the point up to round-off has the key zero; all of them are alike,
  //   - a cell that holds it only within the tolerance has the key by which it misses the point,
  //   - between cells of the same key, the one with the smallest global id wins.
  struct Candidate {
    double key{std::numeric_limits<double>::infinity()};
    GlobalElemId globalId{std::numeric_limits<GlobalElemId>::max()};
  };
  const auto better = [](const Candidate& a, const Candidate& b) {
    return a.key < b.key || (a.key == b.key && a.globalId < b.globalId);
  };

  std::vector<Candidate> best(numPoints);

#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < elements.size(); ++elem) {
    auto planeEquations = std::array<std::array<double, Cell::Dim + 1>, Cell::Dim + 1>();
    auto normLengths = std::array<double, Cell::NumFaces>();
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      VrtxCoords n{};
      VrtxCoords p{};
      MeshTools::pointOnPlane(elements[elem], face, vertices, p);
      MeshTools::normal(elements[elem], face, vertices, n);

      for (std::size_t i = 0; i < Cell::Dim; ++i) {
        planeEquations[i][face] = n[i];
      }
      planeEquations[Cell::Dim][face] = -MeshTools::dot(n, p);
      normLengths[face] = std::sqrt(MeshTools::dot(n, n));
    }

    // The normals are not normalized, their length is twice the area of the face; hence the square
    // root of the largest one is about the edge length of the cell. What lies within a small
    // fraction of it from a face is on that face as far as round-off can tell.
    const double onFace =
        1e-8 * std::sqrt(*std::max_element(normLengths.begin(), normLengths.end()));

    for (std::size_t point = 0; point < numPoints; ++point) {
      // geometric Interpretation (up to numerical errors, hence tolerance parameter):
      // resultFace < 0: The point is inside the face (half-space).
      // resultFace = 0: The point is exactly on the face.
      // resultFace > 0: The point is outside the face.

      // we look for the face, where resultFace is the largest; i.e. the face that will the best
      // "invalidate" our membership in the cell.

      // NOLINTNEXTLINE
      double maxValue = -std::numeric_limits<double>::infinity();
      // NOLINTNEXTLINE
      double maxDistance = -std::numeric_limits<double>::infinity();

#pragma omp simd reduction(max : maxValue, maxDistance)
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {

        double resultFace = 0;
        for (std::size_t dim = 0; dim < Cell::Dim + 1; ++dim) {
          resultFace += planeEquations[dim][face] * points1[point][dim];
        }
        maxValue = std::max(maxValue, resultFace);
        maxDistance = std::max(maxDistance, resultFace / normLengths[face]);
      }

      if (maxValue <= tolerance) {
        // meaning: we're below tolerance to consider the cell
        const Candidate candidate{maxDistance <= onFace ? 0.0 : maxValue, elements[elem].globalId};

#pragma omp critical
        {
          if (better(candidate, best[point])) {
            best[point] = candidate;
            meshIds[point] = static_cast<std::size_t>(elements[elem].localId);
          }
        }
      }
    }
  }

  // now reduce over all ranks for the best fit: first the key, then the global id among the cells
  // that share the best key

  std::vector<double> keys(numPoints);
  for (std::size_t point = 0; point < numPoints; ++point) {
    keys[point] = best[point].key;
  }
  MPI_Allreduce(
      MPI_IN_PLACE, keys.data(), keys.size(), MPI_DOUBLE, MPI_MIN, seissol::Mpi::mpi.comm());

  std::vector<std::uint64_t> globalIds(numPoints);
  for (std::size_t point = 0; point < numPoints; ++point) {
    globalIds[point] = best[point].key == keys[point]
                           ? static_cast<std::uint64_t>(best[point].globalId)
                           : std::numeric_limits<std::uint64_t>::max();
  }
  MPI_Allreduce(MPI_IN_PLACE,
                globalIds.data(),
                globalIds.size(),
                MPI_UINT64_T,
                MPI_MIN,
                seissol::Mpi::mpi.comm());

  std::vector<bool> contained(numPoints);
  for (std::size_t i = 0; i < numPoints; ++i) {
    contained[i] = best[i].key < std::numeric_limits<double>::infinity() &&
                   best[i].key == keys[i] &&
                   static_cast<std::uint64_t>(best[i].globalId) == globalIds[i];
  }
  return contained;
}

} // namespace seissol::initializer

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

CellPlanes::CellPlanes(const std::array<const CoordinateT*, Cell::NumVertices>& vertices) {
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    // as MeshTools::normal and MeshTools::pointOnPlane give them for a cell of a MeshReader
    const auto& first = *vertices[MeshTools::FACE2NODES[face][0]];
    CoordinateT ab{};
    CoordinateT ac{};
    MeshTools::sub(*vertices[MeshTools::FACE2NODES[face][1]], first, ab);
    MeshTools::sub(*vertices[MeshTools::FACE2NODES[face][2]], first, ac);
    CoordinateT n{};
    MeshTools::cross(ab, ac, n);

    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      equations_[i][face] = n[i];
    }
    equations_[Cell::Dim][face] = -MeshTools::dot(n, first);
    normLengths_[face] = std::sqrt(MeshTools::dot(n, n));
  }

  // The normals are not normalized, their length is twice the area of the face; hence the square
  // root of the largest one is about the edge length of the cell. What lies within a small
  // fraction of it from a face is on that face as far as round-off can tell.
  onFace_ = 1e-8 * std::sqrt(*std::ranges::max_element(normLengths_));
}

std::vector<bool> findUniqueMeshIds(const Eigen::Vector3d* points,
                                    const seissol::geometry::MeshReader& mesh,
                                    std::size_t numPoints,
                                    std::size_t* meshIds,
                                    double tolerance) {

  const auto& vertices = mesh.getVertices();
  const auto& elements = mesh.getElements();

  auto points1 = std::vector<CellPlanes::Point>(numPoints);
  for (std::size_t point = 0; point < numPoints; ++point) {
    points1[point] = CellPlanes::homogeneous(points[point]);
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
    std::size_t globalId{std::numeric_limits<std::size_t>::max()};
  };
  const auto better = [](const Candidate& a, const Candidate& b) {
    return a.key < b.key || (a.key == b.key && a.globalId < b.globalId);
  };

  std::vector<Candidate> best(numPoints);

  // an index loop; OpenMP loop constructs take range-based loops only from OpenMP 5.0 on
#pragma omp parallel for schedule(static)
  // NOLINTNEXTLINE(modernize-loop-convert)
  for (std::size_t elem = 0; elem < elements.size(); ++elem) {
    const auto& element = elements[elem];
    const CellPlanes planes({&vertices[element.vertices[0]].coords,
                             &vertices[element.vertices[1]].coords,
                             &vertices[element.vertices[2]].coords,
                             &vertices[element.vertices[3]].coords});

    for (std::size_t point = 0; point < numPoints; ++point) {
      const auto fit = planes.fit(points1[point]);

      if (fit.value <= tolerance) {
        // meaning: we're below tolerance to consider the cell
        const Candidate candidate{.key = planes.key(fit), .globalId = element.globalId};

#pragma omp critical
        {
          if (better(candidate, best[point])) {
            best[point] = candidate;
            meshIds[point] = element.localId;
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
                seissol::Mpi::castToMpiType<std::uint64_t>(),
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

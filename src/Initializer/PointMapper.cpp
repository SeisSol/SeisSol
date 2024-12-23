/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2019, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#include "PointMapper.h"
#include "Parallel/MPI.h"
#include <Geometry/MeshDefinition.h>
#include <Geometry/MeshReader.h>
#include <Geometry/MeshTools.h>
#include <array>
#include <cstring>
#include <mpi.h>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

void findMeshIds(const Eigen::Vector3d* points,
                 const seissol::geometry::MeshReader& mesh,
                 std::size_t numPoints,
                 short* contained,
                 unsigned* meshIds,
                 double tolerance) {
  findMeshIds(
      points, mesh.getVertices(), mesh.getElements(), numPoints, contained, meshIds, tolerance);
}

void findMeshIds(const Eigen::Vector3d* points,
                 const std::vector<Vertex>& vertices,
                 const std::vector<Element>& elements,
                 std::size_t numPoints,
                 short* contained,
                 unsigned* meshIds,
                 double tolerance) {

  memset(contained, 0, numPoints * sizeof(short));

  auto points1 = std::vector<std::array<double, 4>>(numPoints);
  for (std::size_t point = 0; point < numPoints; ++point) {
    points1[point][0] = points[point](0);
    points1[point][1] = points[point](1);
    points1[point][2] = points[point](2);
    points1[point][3] = 1.0;
  }

/// @TODO Could use the code generator for the following
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (std::size_t elem = 0; elem < elements.size(); ++elem) {
    auto planeEquations = std::array<std::array<double, 4>, 4>();
    for (int face = 0; face < 4; ++face) {
      VrtxCoords n{};
      VrtxCoords p{};
      MeshTools::pointOnPlane(elements[elem], face, vertices, p);
      MeshTools::normal(elements[elem], face, vertices, n);

      for (int i = 0; i < 3; ++i) {
        planeEquations[i][face] = n[i];
      }
      planeEquations[3][face] = -MeshTools::dot(n, p);
    }
    for (std::size_t point = 0; point < numPoints; ++point) {
      // NOLINTNEXTLINE
      int notInside = 0;
#ifdef _OPENMP
#pragma omp simd reduction(+ : notInside)
#endif
      for (unsigned face = 0; face < 4; ++face) {
        double resultFace = 0;
        for (unsigned dim = 0; dim < 4; ++dim) {
          resultFace += planeEquations[dim][face] * points1[point][dim];
        }
        notInside += (resultFace > tolerance) ? 1 : 0;
      }
      if (notInside == 0) {
#ifdef _OPENMP
#pragma omp critical
        {
#endif
          /* It might actually happen that a point is found in two tetrahedrons
           * if it lies on the boundary. In this case we arbitrarily assign
           * it to the one with the higher meshId.
           * @todo Check if this is a problem with the numerical scheme. */
          /*if (contained[point] != 0) {
             logError() << "point with id " << point << " was already found in a different
          element!";
          }*/
          const auto localId = static_cast<unsigned>(elements[elem].localId);
          if ((contained[point] == 0) || (meshIds[point] > localId)) {
            contained[point] = 1;
            meshIds[point] = localId;
          }
#ifdef _OPENMP
        }
#endif
      }
    }
  }
}

#ifdef USE_MPI
void cleanDoubles(short* contained, std::size_t numPoints) {
  const auto myrank = seissol::MPI::mpi.rank();
  const auto size = seissol::MPI::mpi.size();

  auto globalContained = std::vector<short>(size * numPoints);
  MPI_Allgather(contained,
                numPoints,
                MPI_SHORT,
                globalContained.data(),
                numPoints,
                MPI_SHORT,
                seissol::MPI::mpi.comm());

  std::size_t cleaned = 0;
  for (std::size_t point = 0; point < numPoints; ++point) {
    if (contained[point] == 1) {
      for (int rank = 0; rank < myrank; ++rank) {
        if (globalContained[rank * numPoints + point] == 1) {
          contained[point] = 0;
          ++cleaned;
          break;
        }
      }
    }
  }

  if (cleaned > 0) {
    logInfo(myrank) << "Cleaned " << cleaned << " double occurring points on rank " << myrank
                    << ".";
  }
}
#endif

} // namespace seissol::initializer

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
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
#include <cstring>
#include <Initializer/MemoryAllocator.h>
#include <utils/logger.h>
#include <Parallel/MPI.h>

void seissol::initializers::findMeshIds(glm::dvec3 const* points, MeshReader const& mesh, unsigned numPoints, short* contained, unsigned* meshIds)
{
  std::vector<Vertex> const& vertices = mesh.getVertices();
  std::vector<Element> const& elements = mesh.getElements();

  memset(contained, 0, numPoints * sizeof(short));

  double (*planeEquations)[4][4] = static_cast<double(*)[4][4]>(seissol::memory::allocate(elements.size() * sizeof(double[4][4]), ALIGNMENT));
  for (unsigned elem = 0; elem < elements.size(); ++elem) {
    for (int face = 0; face < 4; ++face) {
      VrtxCoords n, p;
      MeshTools::pointOnPlane(elements[elem], face, vertices, p);
      MeshTools::normal(elements[elem], face, vertices, n);

      for (unsigned i = 0; i < 3; ++i) {
        planeEquations[elem][i][face] = n[i];
      }
      planeEquations[elem][3][face] = - MeshTools::dot(n, p);
    }
  }

  double (*points1)[4] = new double[numPoints][4];
  for (unsigned point = 0; point < numPoints; ++point) {
    points1[point][0] = points[point].x;
    points1[point][1] = points[point].y;
    points1[point][2] = points[point].z;
    points1[point][3] = 1.0;
  }

/// @TODO Could use the code generator for the following
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
  for (unsigned elem = 0; elem < elements.size(); ++elem) {
#if 0 //defined(__AVX__)
      __m256d zero = _mm256_setzero_pd();
      __m256d planeDims[4];
      for (unsigned i = 0; i < 4; ++i) {
        planeDims[i] = _mm256_load_pd(&planeEquations[elem][i][0]);
      }
#endif
    for (unsigned point = 0; point < numPoints; ++point) {
      int l_notInside = 0;
#if 0 //defined(__AVX__)
      // Not working because <0 => 0 should actually be  <=0 => 0
      /*__m256d result = _mm256_setzero_pd();
      for (unsigned dim = 0; dim < 4; ++dim) {
        result = _mm256_add_pd(result, _mm256_mul_pd(planeDims[dim], _mm256_broadcast_sd(&points1[point][dim])) );
      }
      // >0 => (2^64)-1 ; <0 = 0
      __m256d inside4 = _mm256_cmp_pd(result, zero, _CMP_GE_OQ);
      l_notInside = _mm256_movemask_pd(inside4);*/
#else
      double result[4] = { 0.0, 0.0, 0.0, 0.0 };
      for (unsigned dim = 0; dim < 4; ++dim) {
        for (unsigned face = 0; face < 4; ++face) {
          result[face] += planeEquations[elem][dim][face] * points1[point][dim];
        }
      }
      for (unsigned face = 0; face < 4; ++face) {
        l_notInside += (result[face] > 0.0) ? 1 : 0;
      }
#endif
      if (l_notInside == 0) {
#ifdef _OPENMP
        #pragma omp critical
        {
#endif
          /* It might actually happen that a point is found in two tetrahedrons
           * if it lies on the boundary. In this case we arbitrarily assign
           * it to the one with the higher meshId.
           * @todo Check if this is a problem with the numerical scheme. */
          /*if (contained[point] != 0) {
             logError() << "point with id " << point << " was already found in a different element!";
          }*/
          if (contained[point] == 0 || meshIds[point] > elem) {
            contained[point] = 1;
            meshIds[point] = elem;
          }
#ifdef _OPENMP
        }
#endif
      }
    }
  }

  seissol::memory::free(planeEquations);
  delete[] points1;
}

#ifdef USE_MPI
void seissol::initializers::cleanDoubles(short* contained, unsigned numPoints)
{
  int myrank = seissol::MPI::mpi.rank();
  int size = seissol::MPI::mpi.size();

  short* globalContained = new short[size * numPoints];
  MPI_Allgather(contained, numPoints, MPI_SHORT, globalContained, numPoints, MPI_SHORT, seissol::MPI::mpi.comm());

  unsigned cleaned = 0;
  for (unsigned point = 0; point < numPoints; ++point) {
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
    logInfo(myrank) << "Cleaned " << cleaned << " double occurring points on rank " << myrank << ".";
  }

  delete[] globalContained;
}
#endif

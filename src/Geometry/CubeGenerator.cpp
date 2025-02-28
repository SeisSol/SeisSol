// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CubeGenerator.h"
#include "utils/logger.h"

#include "Parallel/MPI.h"

#include <Geometry/MeshDefinition.h>
#include <Initializer/Parameters/CubeGeneratorParameters.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <iterator>
#include <map>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include "MeshReader.h"

namespace {
using TVertex = std::array<int, 3>;

struct CubeVertex {
  TVertex v;

  bool operator<(const CubeVertex& other) const {
    return (v[0] < other.v[0]) || ((v[0] == other.v[0]) && (v[1] < other.v[1])) ||
           ((v[0] == other.v[0]) && (v[1] == other.v[1]) && (v[2] < other.v[2]));
  }
};

// Index of the vertices of a tetraedra in a cube
// even/odd, index of the tetrahedra, index of vertex, offset of the vertices in x/y/z
const TVertex TetVertices[2][5][4] = {{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                       {{1, 0, 0}, {0, 1, 0}, {1, 1, 1}, {1, 1, 0}},
                                       {{1, 0, 0}, {1, 1, 1}, {0, 0, 1}, {1, 0, 1}},
                                       {{0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {1, 1, 1}},
                                       {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}},
                                      {{{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 0}},
                                       {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}, {1, 0, 0}},
                                       {{0, 0, 0}, {1, 0, 1}, {0, 1, 1}, {0, 0, 1}},
                                       {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
                                       {{0, 0, 0}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}}}};

// neighbor tetrahedra for each face, enumerated per cell (5) and side (4), i.e. 5 * 4 = 20
const int TetSideNeighbors[2][20] = {{3, 3, 3, 0, 1, 3, 0, 2, 2, 2, 2, 1, 0, 1, 3, 1, 3, 0, 0, 2},
                                     {2, 3, 0, 1, 1, 3, 3, 2, 2, 1, 1, 0, 0, 3, 2, 1, 2, 0, 0, 1}};

// orientation of the tetrahedra for each face
const int TetSideOrientations[2][20] = {
    {2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0}};

const char* dim2str(unsigned int dim) {
  switch (dim) {
  case 0:
    return "x";
  case 1:
    return "y";
  case 2:
    return "z";
  default:;
  }

  logError() << "Invalid dimension";
  return "invalid"; // Never reached
}

template <typename A, typename B>
std::pair<B, A> flip_pair(const std::pair<A, B>& p) {
  return std::pair<B, A>(p.second, p.first);
}
} // anonymous namespace

namespace seissol::geometry {

CubeGenerator::CubeGenerator(
    int rank,
    int nProcs,
    const std::string& meshFile,
    const seissol::initializer::parameters::CubeGeneratorParameters& cubeParams)
    : MeshReader(rank), // init base class
      rank(rank), nProcs(nProcs) {
  // get cubeGenerator parameters
  const unsigned int cubeMinX = cubeParams.cubeMinX;
  const unsigned int cubeMaxX = cubeParams.cubeMaxX;
  const unsigned int cubeMinY = cubeParams.cubeMinY;
  const unsigned int cubeMaxY = cubeParams.cubeMaxY;
  const unsigned int cubeMinZ = cubeParams.cubeMinZ;
  const unsigned int cubeMaxZ = cubeParams.cubeMaxZ;
  const unsigned int cubeX = cubeParams.cubeX;
  const unsigned int cubeY = cubeParams.cubeY;
  const unsigned int cubeZ = cubeParams.cubeZ;
  const unsigned int cubePx = cubeParams.cubePx;
  const unsigned int cubePy = cubeParams.cubePy;
  const unsigned int cubePz = cubeParams.cubePz;
  const double cubeScale = cubeParams.cubeS;
  const double cubeScaleX = cubeParams.cubeSx;
  const double cubeScaleY = cubeParams.cubeSy;
  const double cubeScaleZ = cubeParams.cubeSz;
  const double cubeTx = cubeParams.cubeTx;
  const double cubeTy = cubeParams.cubeTy;
  const double cubeTz = cubeParams.cubeTz;

  if (cubePx > 1 && (cubeMinX == 6 || cubeMaxX == 6 || cubeMinY == 6 || cubeMaxY == 6 ||
                     cubeMinZ == 6 || cubeMaxZ == 6)) {
    logWarning()
        << "Atleast one boundary condition is set to 6 (periodic boundary), currently leading "
           "to incorrect results when using more than 1 MPI process";
  }

  // create additional variables necessary for cubeGenerator()
  const std::array<unsigned int, 4> numCubes = {cubeX, cubeY, cubeZ, cubeX * cubeY * cubeZ};
  const std::array<unsigned int, 4> numPartitions = {
      cubePx, cubePy, cubePz, cubePx * cubePy * cubePz};

  // check input arguments
  for (int i = 0; i < 3; i++) {
    if (numCubes[i] < 2) {
      logError() << "Number of cubes in" << dim2str(i) << "dimension must be at least 2";
    }
    if (numCubes[i] % numPartitions[i] != 0) {
      logError() << "Number of cubes in" << dim2str(i) << "dimension can not be distribute to"
                 << numPartitions[i] << "partitions";
    }
    if ((numCubes[i] / numPartitions[i]) % 2 != 0) {
      logError() << "Number of cubes per partition in" << dim2str(i)
                 << "dimension must be a multiple of 2";
    }
    // check if numCubes is multiple of numPartitions, should only fail in numPartitions[0]
    if (numCubes[i] % numPartitions[i] != 0) {
      logError() << "Number of cubes in" << dim2str(i)
                 << "dimenstion must be a multiple of number of threads/processes ="
                 << numPartitions[i];
    }
  }

  // Compute additional sizes
  const std::array<unsigned int, 4> numCubesPerPart = {numCubes[0] / numPartitions[0],
                                                       numCubes[1] / numPartitions[1],
                                                       numCubes[2] / numPartitions[2],
                                                       numCubes[3] / numPartitions[3]};
  const std::array<unsigned long, 4> numElemPerPart = {numCubesPerPart[0] * 5,
                                                       numCubesPerPart[1] * 5,
                                                       numCubesPerPart[2] * 5,
                                                       numCubesPerPart[3] * 5};
  const std::array<unsigned int, 4> numVrtxPerPart = {numCubesPerPart[0] + 1,
                                                      numCubesPerPart[1] + 1,
                                                      numCubesPerPart[2] + 1,
                                                      numCubesPerPart[0] * numCubesPerPart[1] *
                                                          numCubesPerPart[2]};
  const std::array<unsigned int, 3> numBndElements = {2 * numCubesPerPart[1] * numCubesPerPart[2],
                                                      2 * numCubesPerPart[0] * numCubesPerPart[2],
                                                      2 * numCubesPerPart[0] * numCubesPerPart[1]};

  // output file name
  const std::string& fileName = meshFile;

  logInfo() << "Start generating a mesh using the CubeGenerator";
  CubeGenerator::cubeGenerator(numCubes,
                               numPartitions,
                               cubeMinX,
                               cubeMaxX,
                               cubeMinY,
                               cubeMaxY,
                               cubeMinZ,
                               cubeMaxZ,
                               numCubesPerPart,
                               numElemPerPart,
                               numVrtxPerPart,
                               numBndElements,
                               cubeScale,
                               cubeScaleX,
                               cubeScaleY,
                               cubeScaleZ,
                               cubeTx,
                               cubeTy,
                               cubeTz,
                               fileName);
}

void CubeGenerator::cubeGenerator(const std::array<unsigned int, 4> numCubes,
                                  const std::array<unsigned int, 4> numPartitions,
                                  unsigned int boundaryMinx,
                                  unsigned int boundaryMaxx,
                                  unsigned int boundaryMiny,
                                  unsigned int boundaryMaxy,
                                  unsigned int boundaryMinz,
                                  unsigned int boundaryMaxz,
                                  const std::array<unsigned int, 4> numCubesPerPart,
                                  const std::array<unsigned long, 4> numElemPerPart,
                                  const std::array<unsigned int, 4> numVrtxPerPart,
                                  const std::array<unsigned int, 3> numBndElements,
                                  double scale,
                                  double scaleX,
                                  double scaleY,
                                  double scaleZ,
                                  double tx,
                                  double ty,
                                  double tz,
                                  const std::string& meshFile) {

  logInfo() << "Total number of cubes:" << numCubes[0] << 'x' << numCubes[1] << 'x' << numCubes[2]
            << '=' << numCubes[3];
  logInfo() << "Total number of partitions" << numPartitions[0] << 'x' << numPartitions[1] << 'x'
            << numPartitions[2] << '=' << numPartitions[3];
  logInfo() << "Total number of cubes per partition:" << numCubesPerPart[0] << 'x'
            << numCubesPerPart[1] << 'x' << numCubesPerPart[2] << '=' << numCubesPerPart[3];
  logInfo() << "Total number of elements per partition:" << numElemPerPart[0] << 'x'
            << numElemPerPart[1] << 'x' << numElemPerPart[2] << '='
            << numElemPerPart[0] * numElemPerPart[1] * numElemPerPart[2];
  logInfo() << "Using" << omp_get_max_threads() << "threads";

  // Setup MPI Communicator
#ifdef USE_MPI
  MPI_Comm commMaster = MPI_COMM_NULL;
  MPI_Comm_split(seissol::MPI::mpi.comm(), rank % 1 == 0 ? 1 : MPI_UNDEFINED, rank, &commMaster);
#endif // USE_MPI

  size_t bndSize = -1;
  size_t bndElemSize = -1;

  int* sizes = nullptr;
  int maxSize = 0;

  // Get important dimensions
  const size_t partitions = numPartitions[3];

  if (partitions != static_cast<unsigned int>(nProcs)) {
    logError() << "Number of partitions does not match number of MPI ranks.";
  }

  bndSize = 6;
  bndElemSize = *std::max_element(numBndElements.begin(), numBndElements.end());

  // Elements
  sizes = new int[1];
  const int size = numElemPerPart[3];
  sizes[0] = size;
  maxSize = std::max(maxSize, size);
  m_elements.resize(sizes[0]);

  std::vector<CubeVertex> vertices;
  vertices.resize(numElemPerPart[3] * 4);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(3)
#endif // _OPENMP
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        const int odd = (zz + yy + xx) % 2;

        for (unsigned int i = 0; i < 5; i++) {
          for (unsigned int j = 0; j < 4; j++) {
            CubeVertex v{};
            v.v[0] = TetVertices[odd][i][j][0] + xx;
            v.v[1] = TetVertices[odd][i][j][1] + yy;
            v.v[2] = TetVertices[odd][i][j][2] + zz;
            vertices[c] = v;
            c++;
          }
        }
      }
    }
  }

  int* elemVertices = new int[numElemPerPart[3] * 4];
  std::map<CubeVertex, int> vertexMap;

  // Calculate elemVertices
  for (unsigned int i = 0; i < vertices.size(); i++) {
    const auto it = vertexMap.find(vertices[i]);
    if (it != vertexMap.end()) {
      elemVertices[i] = it->second;
    } else {
      const int n = vertexMap.size();
      vertexMap[vertices[i]] = n;
      elemVertices[i] = n;
    }
  }

  int* elemNeighbors = new int[numElemPerPart[3] * 4];
  const int tetNeighbors[2][20] = {
      {-static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5 + 2,
       -static_cast<int>(numCubesPerPart[0]) * 5,
       -4,
       4,
       4,
       -static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5 + 3,
       5,
       static_cast<int>(numCubesPerPart[0]) * 5 + 1,
       4,
       7,
       -static_cast<int>(numCubesPerPart[0]) * 5 + 3,
       static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5 + 1,
       -2,
       static_cast<int>(numCubesPerPart[0]) * 5 + 2,
       4,
       static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5,
       0,
       1,
       2,
       3},
      {-4,
       -static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5 + 3,
       4,
       static_cast<int>(numCubesPerPart[0]) * 5,
       4,
       -static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5 + 2,
       -static_cast<int>(numCubesPerPart[0]) * 5 + 1,
       5,
       4,
       -static_cast<int>(numCubesPerPart[0]) * 5 + 3,
       -3,
       static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5,
       8,
       4,
       static_cast<int>(numCubesPerPart[0]) * 5 + 2,
       static_cast<int>(numCubesPerPart[1] * numCubesPerPart[0]) * 5 + 1,
       0,
       1,
       2,
       3}};

  // Calculate elemNeighbors
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(3)
#endif // _OPENMP
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        const unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        const int odd = (zz + yy + xx) % 2;

        memcpy(&elemNeighbors[c], tetNeighbors[odd], sizeof(int) * 20);
        const int offset = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 5;
        for (int i = 0; i < 20; i++) {
          elemNeighbors[c + i] += offset;
        }

        if (xx == 0) { // first cube in a partition in x dimension
          if (odd != 0) {
            if (boundaryMinx == 6 && numPartitions[0] == 1) {
              elemNeighbors[c] += numCubesPerPart[0] * 5;
              elemNeighbors[c + 10] += numCubesPerPart[0] * 5;
            } else {
              elemNeighbors[c] = numElemPerPart[3];
              elemNeighbors[c + 10] = numElemPerPart[3];
            }
          } else {
            if (boundaryMinx == 6 && numPartitions[0] == 1) {
              elemNeighbors[c + 2] += numCubesPerPart[0] * 5;
              elemNeighbors[c + 12] += numCubesPerPart[0] * 5;
            } else {
              elemNeighbors[c + 2] = numElemPerPart[3];
              elemNeighbors[c + 12] = numElemPerPart[3];
            }
          }
        } else if (xx == numCubesPerPart[0] - 1) { // last cube in a partition in x dimension
          if (odd != 0) {
            if (boundaryMaxx == 6 && numPartitions[0] == 1) {
              elemNeighbors[c + 7] -= numCubesPerPart[0] * 5;
              elemNeighbors[c + 12] -= numCubesPerPart[0] * 5;
            } else {
              elemNeighbors[c + 7] = numElemPerPart[3];
              elemNeighbors[c + 12] = numElemPerPart[3];
            }
          } else {
            if (boundaryMaxx == 6 && numPartitions[0] == 1) {
              elemNeighbors[c + 6] -= numCubesPerPart[0] * 5;
              elemNeighbors[c + 9] -= numCubesPerPart[0] * 5;
            } else {
              elemNeighbors[c + 6] = numElemPerPart[3];
              elemNeighbors[c + 9] = numElemPerPart[3];
            }
          }
        }
        if (yy == 0) { // first cube in a partition in y dimension
          if (odd != 0) {
            if (boundaryMiny == 6 && numPartitions[1] == 1) {
              elemNeighbors[c + 6] += numCubesPerPart[0] * numCubesPerPart[1] * 5;
              elemNeighbors[c + 9] += numCubesPerPart[0] * numCubesPerPart[1] * 5;
            } else {
              elemNeighbors[c + 6] = numElemPerPart[3];
              elemNeighbors[c + 9] = numElemPerPart[3];
            }
          } else {
            if (boundaryMiny == 6 && numPartitions[1] == 1) {
              elemNeighbors[c + 1] += numCubesPerPart[0] * numCubesPerPart[1] * 5;
              elemNeighbors[c + 10] += numCubesPerPart[0] * numCubesPerPart[1] * 5;
            } else {
              elemNeighbors[c + 1] = numElemPerPart[3];
              elemNeighbors[c + 10] = numElemPerPart[3];
            }
          }
        } else if (yy == numCubesPerPart[1] - 1) { // last cube in a partition in y dimension
          if (odd != 0) {
            if (boundaryMaxy == 6 && numPartitions[1] == 1) {
              elemNeighbors[c + 3] -= numCubesPerPart[0] * numCubesPerPart[1] * 5;
              elemNeighbors[c + 14] -= numCubesPerPart[0] * numCubesPerPart[1] * 5;
            } else {
              elemNeighbors[c + 3] = numElemPerPart[3];
              elemNeighbors[c + 14] = numElemPerPart[3];
            }
          } else {
            if (boundaryMaxy == 6 && numPartitions[1] == 1) {
              elemNeighbors[c + 7] -= numCubesPerPart[0] * numCubesPerPart[1] * 5;
              elemNeighbors[c + 13] -= numCubesPerPart[0] * numCubesPerPart[1] * 5;
            } else {
              elemNeighbors[c + 7] = numElemPerPart[3];
              elemNeighbors[c + 13] = numElemPerPart[3];
            }
          }
        }
        if (zz == 0) { // first cube in a partition in z dimension
          if (odd != 0) {
            if (boundaryMinz == 6 && numPartitions[2] == 1) {
              elemNeighbors[c + 1] +=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
              elemNeighbors[c + 5] +=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
            } else {
              elemNeighbors[c + 1] = numElemPerPart[3];
              elemNeighbors[c + 5] = numElemPerPart[3];
            }
          } else {
            if (boundaryMinz == 6 && numPartitions[2] == 1) {
              elemNeighbors[c] += numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
              elemNeighbors[c + 5] +=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
            } else {
              elemNeighbors[c] = numElemPerPart[3];
              elemNeighbors[c + 5] = numElemPerPart[3];
            }
          }
        } else if (zz == numCubesPerPart[2] - 1) { // last cube in a partition in z dimension
          if (odd != 0) {
            if (boundaryMaxz == 6 && numPartitions[2] == 1) {
              elemNeighbors[c + 11] -=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
              elemNeighbors[c + 15] -=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
            } else {
              elemNeighbors[c + 11] = numElemPerPart[3];
              elemNeighbors[c + 15] = numElemPerPart[3];
            }
          } else {
            if (boundaryMaxz == 6 && numPartitions[2] == 1) {
              elemNeighbors[c + 11] -=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
              elemNeighbors[c + 15] -=
                  numCubesPerPart[0] * numCubesPerPart[1] * numCubesPerPart[2] * 5;
            } else {
              elemNeighbors[c + 11] = numElemPerPart[3];
              elemNeighbors[c + 15] = numElemPerPart[3];
            }
          }
        }
      }
    }
  }

  int* elemBoundaries = new int[numElemPerPart[3] * 4];

  // Calculate elemBoundaries
  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      const unsigned int x = rank;
      memset(elemBoundaries, 0, sizeof(int) * numElemPerPart[3] * 4);

      if (x == 0) { // first partition in x dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy) % 2;
            if (odd != 0) {
              elemBoundaries[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] =
                  boundaryMinx;
              elemBoundaries[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 10] =
                  boundaryMinx;
            } else {
              elemBoundaries[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] =
                  boundaryMinx;
              elemBoundaries[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 12] =
                  boundaryMinx;
            }
          }
        }
      }
      if (x == numPartitions[0] - 1) { // last partition in x dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy + 1) % 2;
            if (odd != 0) {
              elemBoundaries[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             7] = boundaryMaxx;
              elemBoundaries[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             12] = boundaryMaxx;
            } else {
              elemBoundaries[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             6] = boundaryMaxx;
              elemBoundaries[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             9] = boundaryMaxx;
            }
          }
        }
      }
      if (y == 0) { // first partition in y dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx) % 2;
            if (odd != 0) {
              elemBoundaries[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 6] =
                  boundaryMiny;
              elemBoundaries[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 9] =
                  boundaryMiny;
            } else {
              elemBoundaries[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 1] =
                  boundaryMiny;
              elemBoundaries[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 10] =
                  boundaryMiny;
            }
          }
        }
      }
      if (y == numPartitions[1] - 1) { // last partition in y dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx + 1) % 2;
            if (odd != 0) {
              elemBoundaries
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   3] = boundaryMaxy;
              elemBoundaries
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   14] = boundaryMaxy;
            } else {
              elemBoundaries
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   7] = boundaryMaxy;
              elemBoundaries
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   13] = boundaryMaxy;
            }
          }
        }
      }
      if (z == 0) { // first partition in z dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (yy + xx) % 2;
            if (odd != 0) {
              elemBoundaries[(yy * numCubesPerPart[0] + xx) * 20 + 1] = boundaryMinz;
              elemBoundaries[(yy * numCubesPerPart[0] + xx) * 20 + 5] = boundaryMinz;
            } else {
              elemBoundaries[(yy * numCubesPerPart[0] + xx) * 20] = boundaryMinz;
              elemBoundaries[(yy * numCubesPerPart[0] + xx) * 20 + 5] = boundaryMinz;
            }
          }
        }
      }
      if (z == numPartitions[2] - 1) { // last partition in z dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            //                                                                    int odd =
            //                                                                    (yy+xx+1) % 2; if
            //                                                                    (odd) {
            elemBoundaries
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 11] = boundaryMaxz;
            elemBoundaries
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 15] = boundaryMaxz;
            //                                                                    } else {
            //                                                                            elemBoundaries[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 11]
            //                                                                            =
            //                                                                            boundaryMaxz;
            //                                                                            elemBoundaries[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 15]
            //                                                                            =
            //                                                                            boundaryMaxz;
            //                                                                    }
          }
        }
      }

      for (int i = 0; i < sizes[0]; i++) {
        // ElemBoundaries is an int array of size 4
        memcpy(m_elements[i].boundaries, &elemBoundaries[i * 4], sizeof(ElemBoundaries));
      }
    }
  }

  int* elemNeighborSides = new int[numElemPerPart[3] * 4];
  int* elemNeighborSidesDef = new int[numElemPerPart[3] * 4];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(3)
#endif // _OPENMP
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        const int odd = (zz + yy + xx) % 2;
        const unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        memcpy(&elemNeighborSidesDef[c], TetSideNeighbors[odd], sizeof(int) * 20);
      }
    }
  }

  // Calculate elemNeighborSides
  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      const unsigned int x = rank;
      memcpy(elemNeighborSides, elemNeighborSidesDef, sizeof(int) * numElemPerPart[3] * 4);

      if (boundaryMinx != 6 && x == 0) { // first partition in x dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy) % 2;
            if (odd != 0) {
              elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] = 0;
              elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 10] = 0;
            } else {
              elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] = 0;
              elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 12] = 0;
            }
          }
        }
      }
      if (boundaryMaxx != 6 && x == numPartitions[0] - 1) { // last partition in x dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy + 1) % 2;
            if (odd != 0) {
              elemNeighborSides[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                7] = 0;
              elemNeighborSides[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                12] = 0;
            } else {
              elemNeighborSides[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                6] = 0;
              elemNeighborSides[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                9] = 0;
            }
          }
        }
      }
      if (boundaryMiny != 6 && y == 0) { // first partition in y dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx) % 2;
            if (odd != 0) {
              elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 6] = 0;
              elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 9] = 0;
            } else {
              elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 1] = 0;
              elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 10] = 0;
            }
          }
        }
      }
      if (boundaryMaxy != 6 && y == numPartitions[1] - 1) { // last partition in y dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx + 1) % 2;
            if (odd != 0) {
              elemNeighborSides
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   3] = 0;
              elemNeighborSides
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   14] = 0;
            } else {
              elemNeighborSides
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   7] = 0;
              elemNeighborSides
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   13] = 0;
            }
          }
        }
      }
      if (boundaryMinz != 6 && z == 0) { // first partition in z dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (yy + xx) % 2;
            if (odd != 0) {
              elemNeighborSides[(yy * numCubesPerPart[0] + xx) * 20 + 1] = 0;
              elemNeighborSides[(yy * numCubesPerPart[0] + xx) * 20 + 5] = 0;
            } else {
              elemNeighborSides[(yy * numCubesPerPart[0] + xx) * 20] = 0;
              elemNeighborSides[(yy * numCubesPerPart[0] + xx) * 20 + 5] = 0;
            }
          }
        }
      }
      if (boundaryMaxz != 6 && z == numPartitions[2] - 1) { // last partition in z dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            //                                                                    int odd =
            //                                                                    (yy+xx+1) % 2; if
            //                                                                    (odd) {
            elemNeighborSides
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 11] = 0;
            elemNeighborSides
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 15] = 0;
            //                                                                    } else {
            //                                                                            elemNeighborSides[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 11]
            //                                                                            = 0;
            //                                                                            elemNeighborSides[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 15]
            //                                                                            = 0;
            //                                                                    }
          }
        }
      }

      for (int i = 0; i < sizes[0]; i++) {
        // ElemNeighborSides is an int array of size 4
        memcpy(m_elements[i].neighborSides, &elemNeighborSides[i * 4], sizeof(ElemNeighborSides));
      }
    }
  }

  delete[] elemNeighborSidesDef;

  int* elemSideOrientations = new int[numElemPerPart[3] * 4];
  int* elemSideOrientationsDef = new int[numElemPerPart[3] * 4];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(3)
#endif // _OPENMP
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        const int odd = (zz + yy + xx) % 2;
        const unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        memcpy(&elemSideOrientationsDef[c], TetSideOrientations[odd], sizeof(int) * 20);
      }
    }
  }

  // Calculate elemSideOrientations
  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      const unsigned int x = rank;

      memcpy(elemSideOrientations, elemSideOrientationsDef, sizeof(int) * numElemPerPart[3] * 4);

      if (boundaryMinx != 6 && x == 0) { // first partition in x dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy) % 2;
            if (odd != 0) {
              elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] = 0;
              elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 10] =
                  0;
            } else {
              elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] =
                  0;
              elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 12] =
                  0;
            }
          }
        }
      }
      if (boundaryMaxx != 6 && x == numPartitions[0] - 1) { // last partition in x dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy + 1) % 2;
            if (odd != 0) {
              elemSideOrientations[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                    numCubesPerPart[0] - 1) *
                                       20 +
                                   7] = 0;
              elemSideOrientations[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                    numCubesPerPart[0] - 1) *
                                       20 +
                                   12] = 0;
            } else {
              elemSideOrientations[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                    numCubesPerPart[0] - 1) *
                                       20 +
                                   6] = 0;
              elemSideOrientations[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                    numCubesPerPart[0] - 1) *
                                       20 +
                                   9] = 0;
            }
          }
        }
      }
      // There are zero anyway
      //                              if (boundaryMiny != 6 && y == 0) { // first partition in y
      //                              dimension
      //                                      #ifdef _OPENMP
      //                                      #pragma omp parallel for schedule(static) collapse(2)
      //                                      #endig // _OPENMP
      //                                      for (unsigned int zz = 0; zz < numCubesPerPart[2];
      //                                      zz++) {
      //                                              for (unsigned int xx = 0; xx <
      //                                              numCubesPerPart[0]; xx++) {
      //                                                      int odd = (zz+xx) % 2;
      //                                                      if (odd) {
      //                                                              elemSideOrientations[(zz*numCubesPerPart[1]*numCubesPerPart[0]+xx)
      //                                                              * 20 + 6] = 0;
      //                                                              elemSideOrientations[(zz*numCubesPerPart[1]*numCubesPerPart[0]+xx)
      //                                                              * 20 + 9] = 0;
      //                                                      } else {
      //                                                              elemSideOrientations[(zz*numCubesPerPart[1]*numCubesPerPart[0]+xx)
      //                                                              * 20 + 1] = 0;
      //                                                              elemSideOrientations[(zz*numCubesPerPart[1]*numCubesPerPart[0]+xx)
      //                                                              * 20 + 10] = 0;
      //                                                      }
      //                                              }
      //                                      }
      //                              }
      //                              if (boundaryMaxy != 6 && y == numPartitions[1]-1) { //
      //                              last partition in y dimension
      //                                      #ifdef _OPENMP
      //                                      #pragma omp parallel for schedule(static) collapse(2)
      //                                      #endig // _OPENMP
      //                                      for (unsigned int zz = 0; zz < numCubesPerPart[2];
      //                                      zz++) {
      //                                              for (unsigned int xx = 0; xx <
      //                                              numCubesPerPart[0]; xx++) {
      //                                                      int odd = (zz+xx+1) % 2;
      //                                                      if (odd) {
      //                                                              elemSideOrientations[((zz*numCubesPerPart[1]+numCubesPerPart[1]-1)*numCubesPerPart[0]+xx)
      //                                                              * 20 + 3] = 0;
      //                                                              elemSideOrientations[((zz*numCubesPerPart[1]+numCubesPerPart[1]-1)*numCubesPerPart[0]+xx)
      //                                                              * 20 + 14] = 0;
      //                                                      } else {
      //                                                              elemSideOrientations[((zz*numCubesPerPart[1]+numCubesPerPart[1]-1)*numCubesPerPart[0]+xx)
      //                                                              * 20 + 7] = 0;
      //                                                              elemSideOrientations[((zz*numCubesPerPart[1]+numCubesPerPart[1]-1)*numCubesPerPart[0]+xx)
      //                                                              * 20 + 13] = 0;
      //                                                      }
      //                                              }
      //                                      }
      //                              }
      if (boundaryMinz != 6 && z == 0) { // first partition in z dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (yy + xx) % 2;
            if (odd != 0) {
              elemSideOrientations[(yy * numCubesPerPart[0] + xx) * 20 + 1] = 0;
              elemSideOrientations[(yy * numCubesPerPart[0] + xx) * 20 + 5] = 0;
            } else {
              elemSideOrientations[(yy * numCubesPerPart[0] + xx) * 20] = 0;
              elemSideOrientations[(yy * numCubesPerPart[0] + xx) * 20 + 5] = 0;
            }
          }
        }
      }
      if (boundaryMaxz != 6 && z == numPartitions[2] - 1) { // last partition in z dimension
#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            //                                                                    int odd =
            //                                                                    (yy+xx+1) % 2; if
            //                                                                    (odd) {
            elemSideOrientations
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 11] = 0;
            elemSideOrientations
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 15] = 0;
            //                                                                    } else {
            //                                                                            elemSideOrientations[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 11]
            //                                                                            = 0;
            //                                                                            elemSideOrientations[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 15]
            //                                                                            = 0;
            //                                                                    }
          }
        }
      }

      for (int i = 0; i < sizes[0]; i++) {
        // ElemSideOrientations is an int array of size 4
        memcpy(m_elements[i].sideOrientations,
               &elemSideOrientations[i * 4],
               sizeof(ElemSideOrientations));
      }
    }
  }

  delete[] elemSideOrientationsDef;

  int* elemNeighborRanks = new int[numElemPerPart[3] * 4];

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      const unsigned int x = rank;
      const int myrank = (z * numPartitions[1] + y) * numPartitions[0] + x;

      std::fill(elemNeighborRanks, elemNeighborRanks + numElemPerPart[3] * 4, myrank);

      if ((boundaryMinx == 6 && numPartitions[0] > 1) || x != 0) { // first partition in x dimension
        const int rank = (z * numPartitions[1] + y) * numPartitions[0] +
                         (x - 1 + numPartitions[0]) % numPartitions[0];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy) % 2;
            if (odd != 0) {
              elemNeighborRanks[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] = rank;
              elemNeighborRanks[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 10] =
                  rank;
            } else {
              elemNeighborRanks[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] =
                  rank;
              elemNeighborRanks[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 12] =
                  rank;
            }
          }
        }
      }
      if ((boundaryMaxx == 6 && numPartitions[0] > 1) ||
          x != numPartitions[0] - 1) { // last partition in x dimension
        const int rank = (z * numPartitions[1] + y) * numPartitions[0] + (x + 1) % numPartitions[0];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy + 1) % 2;
            if (odd != 0) {
              elemNeighborRanks[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                7] = rank;
              elemNeighborRanks[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                12] = rank;
            } else {
              elemNeighborRanks[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                6] = rank;
              elemNeighborRanks[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                                 numCubesPerPart[0] - 1) *
                                    20 +
                                9] = rank;
            }
          }
        }
      }
      if ((boundaryMiny == 6 && numPartitions[1] > 1) || y != 0) { // first partition in y dimension
        const int rank = (z * numPartitions[1] + (y - 1 + numPartitions[1]) % numPartitions[1]) *
                             numPartitions[0] +
                         x;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx) % 2;
            if (odd != 0) {
              elemNeighborRanks[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 6] =
                  rank;
              elemNeighborRanks[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 9] =
                  rank;
            } else {
              elemNeighborRanks[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 1] =
                  rank;
              elemNeighborRanks[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 10] =
                  rank;
            }
          }
        }
      }
      if ((boundaryMaxy == 6 && numPartitions[1] > 1) ||
          y != numPartitions[1] - 1) { // last partition in y dimension
        const int rank = (z * numPartitions[1] + (y + 1) % numPartitions[1]) * numPartitions[0] + x;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx + 1) % 2;
            if (odd != 0) {
              elemNeighborRanks
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   3] = rank;
              elemNeighborRanks
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   14] = rank;
            } else {
              elemNeighborRanks
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   7] = rank;
              elemNeighborRanks
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   13] = rank;
            }
          }
        }
      }
      if ((boundaryMinz == 6 && numPartitions[2] > 1) || z != 0) { // first partition in z dimension
        const int rank = (((z - 1 + numPartitions[2]) % numPartitions[2]) * numPartitions[1] + y) *
                             numPartitions[0] +
                         x;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (yy + xx) % 2;
            if (odd != 0) {
              elemNeighborRanks[(yy * numCubesPerPart[0] + xx) * 20 + 1] = rank;
              elemNeighborRanks[(yy * numCubesPerPart[0] + xx) * 20 + 5] = rank;
            } else {
              elemNeighborRanks[(yy * numCubesPerPart[0] + xx) * 20] = rank;
              elemNeighborRanks[(yy * numCubesPerPart[0] + xx) * 20 + 5] = rank;
            }
          }
        }
      }
      if ((boundaryMaxz == 6 && numPartitions[2] > 1) ||
          z != numPartitions[2] - 1) { // last partition in z dimension
        const int rank =
            (((z + 1) % numPartitions[2]) * numPartitions[1] + y) * numPartitions[0] + x;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) collapse(2)
#endif // _OPENMP
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            //                                                                    int odd =
            //                                                                    (yy+xx+1) % 2; if
            //                                                                    (odd) {
            elemNeighborRanks
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 11] = rank;
            elemNeighborRanks
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 15] = rank;
            //                                                                    } else {
            //                                                                            elemNeighborRanks[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 11]
            //                                                                            = rank;
            //                                                                            elemNeighborRanks[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
            //                                                                            * 20 + 15]
            //                                                                            = rank;
            //                                                                    }
          }
        }
      }

      for (int i = 0; i < sizes[0]; i++) {
        // ElemNeighborRanks is an int array of size 4
        memcpy(m_elements[i].neighborRanks, &elemNeighborRanks[i * 4], sizeof(ElemNeighborRanks));
      }
    }
  }

  int* elemMPIIndices = new int[numElemPerPart[3] * 4];
  int* bndLocalIds = new int[*std::max_element(numBndElements.begin(), numBndElements.end())];

  auto* bndSizePtr = new size_t[numPartitions[3]];
  int* bndElemSizePtr = new int[numPartitions[3] * bndSize];
  int* bndElemRankPtr = new int[numPartitions[3] * bndSize];
  int* bndElemLocalIdsPtr = new int[numPartitions[3] * bndSize * bndElemSize];
  int* elemMPIIndicesPtr = new int[numPartitions[3] * numElemPerPart[3] * 4];

  const int bndSizeGlobal = bndSize;

  // calculate bndElem variables, bndLocalIds and elemMPIIndices
  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      const unsigned int x = rank;
      memset(elemMPIIndices, 0, sizeof(int) * numElemPerPart[3] * 4);

      unsigned int bndSize = 0;

      if ((boundaryMinz == 6 && numPartitions[2] > 1) || z != 0) { // first partition in z dimension
        int nextMPIIndex = 0;
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (yy + xx) % 2;
            if (odd != 0) {
              bndLocalIds[nextMPIIndex] = (yy * numCubesPerPart[0] + xx) * 5 + 1;
              elemMPIIndices[(yy * numCubesPerPart[0] + xx) * 20 + 5] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] = (yy * numCubesPerPart[0] + xx) * 5;
              elemMPIIndices[(yy * numCubesPerPart[0] + xx) * 20 + 1] = nextMPIIndex++;
            } else {
              bndLocalIds[nextMPIIndex] = (yy * numCubesPerPart[0] + xx) * 5;
              elemMPIIndices[(yy * numCubesPerPart[0] + xx) * 20] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] = (yy * numCubesPerPart[0] + xx) * 5 + 1;
              elemMPIIndices[(yy * numCubesPerPart[0] + xx) * 20 + 5] = nextMPIIndex++;
            }
          }
        }

        const size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0U};
        const size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
        const int rank = (((z - 1 + numPartitions[2]) % numPartitions[2]) * numPartitions[1] + y) *
                             numPartitions[0] +
                         x;

        bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
        bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
        memcpy(&bndElemLocalIdsPtr[(start[0] * bndSizeGlobal + start[1]) * bndElemSize + start[2]],
               bndLocalIds,
               sizeof(int) * count[2]);

        bndSize++;
      }
      if ((boundaryMiny == 6 && numPartitions[1] > 1) || y != 0) { // first partition in y dimension
        int nextMPIIndex = 0;
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx) % 2;
            if (odd != 0) {
              bndLocalIds[nextMPIIndex] =
                  (zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 5 + 1;
              elemMPIIndices[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 6] =
                  nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  (zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 5 + 2;
              elemMPIIndices[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 9] =
                  nextMPIIndex++;
            } else {
              bndLocalIds[nextMPIIndex] = (zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 5;
              elemMPIIndices[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 1] =
                  nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  (zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 5 + 2;
              elemMPIIndices[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 10] =
                  nextMPIIndex++;
            }
          }
        }

        const size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0U};
        const size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
        const int rank = (z * numPartitions[1] + (y - 1 + numPartitions[1]) % numPartitions[1]) *
                             numPartitions[0] +
                         x;

        bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
        bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
        memcpy(&bndElemLocalIdsPtr[(start[0] * bndSizeGlobal + start[1]) * bndElemSize + start[2]],
               bndLocalIds,
               sizeof(int) * count[2]);

        bndSize++;
      }
      if ((boundaryMinx == 6 && numPartitions[0] > 1) || x != 0) { // first partition in x dimension
        int nextMPIIndex = 0;
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy) % 2;
            if (odd != 0) {
              bndLocalIds[nextMPIIndex] = (zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 5;
              elemMPIIndices[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] =
                  nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  (zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 5 + 2;
              elemMPIIndices[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 10] =
                  nextMPIIndex++;
            } else {
              bndLocalIds[nextMPIIndex] = (zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 5;
              elemMPIIndices[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] =
                  nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  (zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 5 + 3;
              elemMPIIndices[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 12] =
                  nextMPIIndex++;
            }
          }
        }

        const size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0U};
        const size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
        const int rank = (z * numPartitions[1] + y) * numPartitions[0] +
                         (x - 1 + numPartitions[0]) % numPartitions[0];

        bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
        bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
        memcpy(&bndElemLocalIdsPtr[(start[0] * bndSizeGlobal + start[1]) * bndElemSize + start[2]],
               bndLocalIds,
               sizeof(int) * count[2]);

        bndSize++;
      }
      if ((boundaryMaxx == 6 && numPartitions[0] > 1) ||
          x != numPartitions[0] - 1) { // last partition in x dimension
        int nextMPIIndex = 0;
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            const int odd = (zz + yy + 1) % 2;
            if (odd != 0) {
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + numCubesPerPart[0] - 1) *
                      5 +
                  1;
              elemMPIIndices[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             7] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + numCubesPerPart[0] - 1) *
                      5 +
                  3;
              elemMPIIndices[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             12] = nextMPIIndex++;
            } else {
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + numCubesPerPart[0] - 1) *
                      5 +
                  1;
              elemMPIIndices[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             6] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + numCubesPerPart[0] - 1) *
                      5 +
                  2;
              elemMPIIndices[((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] +
                              numCubesPerPart[0] - 1) *
                                 20 +
                             9] = nextMPIIndex++;
            }
          }
        }

        const size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0};
        const size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
        int rank = (z * numPartitions[1] + y) * numPartitions[0] + (x + 1) % numPartitions[0];
        rank = (rank + numPartitions[3]) % numPartitions[3];

        bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
        bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
        memcpy(&bndElemLocalIdsPtr[(start[0] * bndSizeGlobal + start[1]) * bndElemSize + start[2]],
               bndLocalIds,
               sizeof(int) * count[2]);

        bndSize++;
      }
      if ((boundaryMaxy == 6 && numPartitions[1] > 1) ||
          y != numPartitions[1] - 1) { // last partition in y dimension
        int nextMPIIndex = 0;
        for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            const int odd = (zz + xx + 1) % 2;
            if (odd != 0) {
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                  5;
              elemMPIIndices
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   3] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                      5 +
                  3;
              elemMPIIndices
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   14] = nextMPIIndex++;
            } else {
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                      5 +
                  1;
              elemMPIIndices
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   7] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                      5 +
                  3;
              elemMPIIndices
                  [((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                       20 +
                   13] = nextMPIIndex++;
            }
          }
        }

        const size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0};
        const size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
        int rank = (z * numPartitions[1] + (y + 1) % numPartitions[1]) * numPartitions[0] + x;
        rank = (rank + numPartitions[3]) % numPartitions[3];

        bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
        bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
        memcpy(&bndElemLocalIdsPtr[(start[0] * bndSizeGlobal + start[1]) * bndElemSize + start[2]],
               bndLocalIds,
               sizeof(int) * count[2]);

        bndSize++;
      }
      if ((boundaryMaxz == 6 && numPartitions[2] > 1) ||
          z != numPartitions[2] - 1) { // last partition in z dimension
        int nextMPIIndex = 0;
        for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
          for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
            bndLocalIds[nextMPIIndex] =
                (((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                    5 +
                2;
            elemMPIIndices
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 11] = nextMPIIndex++;
            bndLocalIds[nextMPIIndex] =
                (((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                    5 +
                3;
            elemMPIIndices
                [(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                     20 +
                 15] = nextMPIIndex++;
          }
        }

        const size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0U};
        const size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
        int rank = (((z + 1) % numPartitions[2]) * numPartitions[1] + y) * numPartitions[0] + x;
        rank = (rank + numPartitions[3]) % numPartitions[3];

        bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
        bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
        memcpy(&bndElemLocalIdsPtr[(start[0] * bndSizeGlobal + start[1]) * bndElemSize + start[2]],
               bndLocalIds,
               sizeof(int) * count[2]);

        bndSize++;
      }

      for (int i = 0; i < sizes[0]; i++) {
        // ElemMPIIndices is an int array of size 4
        memcpy(m_elements[i].mpiIndices, &elemMPIIndices[i * 4], sizeof(ElemMPIIndices));
      }

      bndSizePtr[(z * numPartitions[1] + y) * numPartitions[0] + x] = bndSize;
    }
  }

  // Set material zone to 1
  int* elemGroup = new int[numElemPerPart[3]];
  std::fill(elemGroup, elemGroup + numElemPerPart[3], 1);

  // copy the remaining Elem variables to m_elements
  auto* elemVerticesCast = reinterpret_cast<ElemVertices*>(elemVertices);
  auto* elemNeighborsCast = reinterpret_cast<ElemNeighbors*>(elemNeighbors);

  for (int i = 0; i < sizes[0]; i++) {
    m_elements[i].localId = i;

    memcpy(m_elements[i].vertices, &elemVerticesCast[i], sizeof(ElemVertices));
    memcpy(m_elements[i].neighbors, &elemNeighborsCast[i], sizeof(ElemNeighbors));
    m_elements[i].group = elemGroup[i];
  }

  delete[] elemVerticesCast;
  delete[] elemNeighborsCast;
  delete[] elemNeighborSides;
  delete[] elemSideOrientations;
  delete[] elemBoundaries;
  delete[] elemNeighborRanks;
  delete[] elemMPIIndices;
  delete[] elemGroup;

  // Vertices
  std::map<int, CubeVertex> uniqueVertices;
  transform(vertexMap.begin(),
            vertexMap.end(),
            inserter(uniqueVertices, uniqueVertices.begin()),
            flip_pair<CubeVertex, int>);

  m_vertices.resize(uniqueVertices.size());

  const double halfWidthX = scaleX / 2.0;
  const double halfWidthY = scaleY / 2.0;
  const double halfWidthZ = scaleZ / 2.0;

  // Calculate vrtxCoords
  auto* vrtxCoords = new double[uniqueVertices.size() * 3];

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      const unsigned int x = rank;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif // _OPENMP
      for (unsigned int i = 0; i < uniqueVertices.size(); i++) {
        vrtxCoords[i * 3] =
            static_cast<double>(uniqueVertices.at(i).v[0] + x * numCubesPerPart[0]) /
                static_cast<double>(numCubes[0]) * scaleX -
            halfWidthX + tx;
        vrtxCoords[i * 3 + 1] =
            static_cast<double>(uniqueVertices.at(i).v[1] + y * numCubesPerPart[1]) /
                static_cast<double>(numCubes[1]) * scaleY -
            halfWidthY + ty;
        vrtxCoords[i * 3 + 2] =
            static_cast<double>(uniqueVertices.at(i).v[2] + z * numCubesPerPart[2]) /
                static_cast<double>(numCubes[2]) * scaleZ -
            halfWidthZ + tz;
      }
    }
  }

  // Copy buffers to vertices
  for (std::size_t i = 0; i < uniqueVertices.size(); i++) {
    // VrtxCoord is defined as an int array of size 3
    memcpy(m_vertices[i].coords, &vrtxCoords[i * 3], sizeof(VrtxCoords));
  }

  delete[] vrtxCoords;

  // Adjust sizes for the upcoming determination of neighbors
  sizes[0] = bndSizePtr[static_cast<size_t>(rank)];

  // Get maximum number of neighbors (required to get collective MPI-IO right)
  const int maxNeighbors = bndSize;
  // MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, seissol::MPI::mpi.comm());
  int* bndElemLocalIds = new int[bndElemSize];

  //        SCOREP_USER_REGION_DEFINE( r_read_boundaries );
  //        SCOREP_USER_REGION_BEGIN( r_read_boundaries, "read_boundaries",
  //        SCOREP_USER_REGION_TYPE_COMMON );

  size_t bndStart[3] = {0, 0, 0};
  for (int i = 0; i < maxNeighbors; i++) {
    bndStart[0] = static_cast<size_t>(rank);
    bndStart[1] = static_cast<size_t>(i);

    // Get neighbor rank
    const int bndRank = bndElemRankPtr[bndStart[0] * bndSize + bndStart[1]];

    // Get size of this boundary
    const int elemSize = bndElemSizePtr[bndStart[0] * bndSize + bndStart[1]];

    const size_t bndCount[3] = {1, 1, bndElemSize};
    memcpy(bndElemLocalIds,
           &bndElemLocalIdsPtr[(bndStart[0] * bndSize + bndStart[1]) * bndElemSize],
           sizeof(int) * bndCount[2]);

    if (i < sizes[0]) {
      addMPINeighbor(i, bndRank, elemSize, bndElemLocalIds);
    }
  }

  delete[] bndLocalIds;
  delete[] bndElemLocalIds;

  //        SCOREP_USER_REGION_END( r_read_boundaries )

  delete[] sizes;

  delete[] bndSizePtr;
  delete[] bndElemSizePtr;
  delete[] bndElemRankPtr;
  delete[] bndElemLocalIdsPtr;
  delete[] elemMPIIndicesPtr;

  // Close MPI communicator
#ifdef USE_MPI
  MPI_Comm_free(&commMaster);
#endif // USE_MPI

  // Recompute additional information
  findElementsPerVertex();

  logInfo() << "Finished";
}

void CubeGenerator::findElementsPerVertex() {
  for (auto i = m_elements.begin(); i != m_elements.end(); i++) {
    for (int j = 0; j < 4; j++) {
      assert(i->vertices[j] < static_cast<int>(m_vertices.size()));
      // push back the localIds for each element of a vertex
      m_vertices[i->vertices[j]].elements.push_back(i->localId);
    }
  }
}

void CubeGenerator::addMPINeighbor(int localID,
                                   int bndRank,
                                   int elemSize,
                                   const int* bndElemLocalIds) {

  MPINeighbor neighbor;
  neighbor.localID = localID;

  neighbor.elements.resize(elemSize);

  for (int i = 0; i < elemSize; i++) {
    neighbor.elements[i].localElement = bndElemLocalIds[i];
  }

  m_MPINeighbors[bndRank] = neighbor;
}

bool CubeGenerator::inlineTimestepCompute() const { return false; }

bool CubeGenerator::inlineClusterCompute() const { return false; }

} // namespace seissol::geometry

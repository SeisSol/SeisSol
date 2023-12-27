#include "CubeGenerator.h"

#ifdef USE_NETCDF
#include <algorithm>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <netcdf.h>
#include <omp.h>

#include "MeshReader.h"
#include "utils/args.h"
#include "utils/logger.h"


namespace {
typedef int t_vertex[3];

struct CubeVertex {
  t_vertex v;

  bool operator<(const CubeVertex& other) const {
    return (v[0] < other.v[0]) || ((v[0] == other.v[0]) && (v[1] < other.v[1])) ||
           ((v[0] == other.v[0]) && (v[1] == other.v[1]) && (v[2] < other.v[2]));
  }
};

// Index of the vertices of a tetraedra in a cube
// even/odd, index of the tetrahedra, index of vertex, offset of the vertices in x/y/z
static const t_vertex TET_VERTICES[2][5][4] = {{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                                {{1, 0, 0}, {0, 1, 0}, {1, 1, 1}, {1, 1, 0}},
                                                {{1, 0, 0}, {1, 1, 1}, {0, 0, 1}, {1, 0, 1}},
                                                {{0, 1, 0}, {0, 1, 1}, {0, 0, 1}, {1, 1, 1}},
                                                {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}},
                                               {{{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 0}},
                                                {{0, 0, 0}, {1, 1, 0}, {1, 0, 1}, {1, 0, 0}},
                                                {{0, 0, 0}, {1, 0, 1}, {0, 1, 1}, {0, 0, 1}},
                                                {{1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
                                                {{0, 0, 0}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}}}};

static const int TET_SIDE_NEIGHBORS[2][5 * 4] = {
    {3, 3, 3, 0, 1, 3, 0, 2, 2, 2, 2, 1, 0, 1, 3, 1, 3, 0, 0, 2},
    {2, 3, 0, 1, 1, 3, 3, 2, 2, 1, 1, 0, 0, 3, 2, 1, 2, 0, 0, 1}};

static const int TET_SIDE_ORIENTATIONS[2][5 * 4] = {
    {2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0}};
} // anonymous namespace

static inline void loadBar(int x, int n, int r = 100, int w = 50, int rank = 0) {
  std::ostringstream bar;
  r = std::min(r, n);

  // Only update r times.
  if (x % (n / r) != 0)
    return;

  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x / (float)n;
  int c = ratio * w;

  // Show the percentage complete.
  bar << std::setw(3) << (int)(ratio * 100) << "% [";

  // Show the load bar.
  for (int x = 0; x < c; x++)
    bar << '=';

  for (int x = c; x < w; x++)
    bar << ' ';

  bar << ']';

  logInfo(rank) << bar.str();
}

static const char* dim2str(unsigned int dim) {
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
static std::pair<B, A> flip_pair(const std::pair<A, B>& p) {
  return std::pair<B, A>(p.second, p.first);
}

static void checkNcError(int error) {
  if (error != NC_NOERR)
    logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

seissol::geometry::CubeGenerator::CubeGenerator(
    int rank,
    int nProcs,
    const std::string& meshFile,
    const seissol::initializers::parameters::CubeGeneratorParameters& cubeParams)
    : seissol::geometry::MeshReader(rank) {
  // get cubeGenerator parameters
  unsigned int cubeMinX = cubeParams.cubeMinX;
  unsigned int cubeMaxX = cubeParams.cubeMaxX;
  unsigned int cubeMinY = cubeParams.cubeMinY;
  unsigned int cubeMaxY = cubeParams.cubeMaxY;
  unsigned int cubeMinZ = cubeParams.cubeMinZ;
  unsigned int cubeMaxZ = cubeParams.cubeMaxZ;
  unsigned int cubeX = cubeParams.cubeX;
  unsigned int cubeY = cubeParams.cubeY;
  unsigned int cubeZ = cubeParams.cubeZ;
  unsigned int cubePx = cubeParams.cubePx;
  unsigned int cubePy = cubeParams.cubePy;
  unsigned int cubePz = cubeParams.cubePz;
  double cubeScale = cubeParams.cubeS;
  double cubeScaleX = cubeParams.cubeSx;
  double cubeScaleY = cubeParams.cubeSy;
  double cubeScaleZ = cubeParams.cubeSz;
  double cubeTx = cubeParams.cubeTx;
  double cubeTy = cubeParams.cubeTy;
  double cubeTz = cubeParams.cubeTz;

  // create additional variables necessary for cubeGenerator()
  unsigned int numCubes[4];
  unsigned int numPartitions[4];

  // check input arguments
  numCubes[0] = cubeX;
  numCubes[1] = cubeY;
  numCubes[2] = cubeZ;

  numPartitions[0] = cubePx;
  numPartitions[1] = cubePy;
  numPartitions[2] = cubePz;

  for (int i = 0; i < 3; i++) {
    if (numCubes[i] < 2)
      logError() << "Number of cubes in" << dim2str(i) << "dimension must be at least 2";
    if (numCubes[i] % numPartitions[i] != 0)
      logError() << "Number of cubes in" << dim2str(i) << "dimension can not be distribute to"
                 << numPartitions[i] << "partitions";
    if ((numCubes[i] / numPartitions[i]) % 2 != 0)
      logError() << "Number of cubes per partition in" << dim2str(i)
                 << "dimension must be a multiple of 2";
    // check if numCubes is multiple of numPartitions, can only fail in numPartitions[0]
    if (numCubes[i] % numPartitions[i] != 0)
      logError() << "Number of cubes in" << dim2str(i)
                 << "dimenstion must be a multiple of number of threads/processes ="
                 << numPartitions[i];
  }

  // Compute additional sizes
  numCubes[3] = numCubes[0] * numCubes[1] * numCubes[2];
  numPartitions[3] = numPartitions[0] * numPartitions[1] * numPartitions[2];

  unsigned int numCubesPerPart[4];
  unsigned long numElemPerPart[4];
  for (int i = 0; i < 4; i++) {
    numCubesPerPart[i] = numCubes[i] / numPartitions[i];
    numElemPerPart[i] = numCubesPerPart[i] * 5;
  }

  unsigned int numVrtxPerPart[4];
  for (int i = 0; i < 3; i++)
    numVrtxPerPart[i] = numCubesPerPart[i] + 1;
  numVrtxPerPart[3] = numVrtxPerPart[0] * numVrtxPerPart[1] * numVrtxPerPart[2];

  unsigned int numBndElements[3];
  numBndElements[0] = 2 * numCubesPerPart[1] * numCubesPerPart[2];
  numBndElements[1] = 2 * numCubesPerPart[0] * numCubesPerPart[2];
  numBndElements[2] = 2 * numCubesPerPart[0] * numCubesPerPart[1];

  // output file name
  std::string fileName = meshFile;

  logInfo(rank) << "Start generating a mesh using the CubeGenerator";
  seissol::geometry::CubeGenerator::cubeGenerator(numCubes,
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
                                                  fileName.c_str());
}

void seissol::geometry::CubeGenerator::cubeGenerator(unsigned int numCubes[4],
                                                     unsigned int numPartitions[4],
                                                     unsigned int boundaryMinx,
                                                     unsigned int boundaryMaxx,
                                                     unsigned int boundaryMiny,
                                                     unsigned int boundaryMaxy,
                                                     unsigned int boundaryMinz,
                                                     unsigned int boundaryMaxz,
                                                     unsigned int numCubesPerPart[4],
                                                     unsigned long numElemPerPart[4],
                                                     unsigned int numVrtxPerPart[4],
                                                     unsigned int numBndElements[3],
                                                     double scale,
                                                     double scaleX,
                                                     double scaleY,
                                                     double scaleZ,
                                                     double tx,
                                                     double ty,
                                                     double tz,
                                                     const std::string& output) {

  logInfo() << "Total number of cubes:" << numCubes[0] << 'x' << numCubes[1] << 'x' << numCubes[2]
            << '=' << numCubes[3];
  logInfo() << "Total number of partitions" << numPartitions[0] << 'x' << numPartitions[1] << 'x'
            << numPartitions[2] << '=' << numPartitions[3];
  logInfo() << "Total number of cubes per partition:" << numCubesPerPart[0] << 'x'
            << numCubesPerPart[1] << 'x' << numCubesPerPart[2] << '=' << numCubesPerPart[3];
  logInfo() << "Total number of elements per partition:" << numElemPerPart[0] << 'x'
            << numElemPerPart[1] << 'x' << numElemPerPart[2] << '=' << numElemPerPart[3];
  logInfo() << "Using" << omp_get_max_threads() << "threads";

  int netcdf_writes = 2 + numPartitions[3] * 8;

  // Create the netcdf file
  int ncFile;

  checkNcError(nc_create(output.c_str(), NC_NETCDF4, &ncFile));

  int ncDimDimension;
  nc_def_dim(ncFile, "dimension", 3, &ncDimDimension);

  int ncDimPart;
  nc_def_dim(ncFile, "partitions", numPartitions[3], &ncDimPart);

  int ncDimElem, ncDimElemSides, ncDimElemVertices;
  nc_def_dim(ncFile, "elements", numElemPerPart[3], &ncDimElem);
  nc_def_dim(ncFile, "element_sides", 4, &ncDimElemSides);
  nc_def_dim(ncFile, "element_vertices_dim", 4, &ncDimElemVertices);

  int ncDimVrtx;
  nc_def_dim(ncFile, "vertices", numVrtxPerPart[3], &ncDimVrtx);

  int ncDimBnd, ncDimBndElem;
  nc_def_dim(ncFile, "boundaries", 6, &ncDimBnd);
  nc_def_dim(ncFile,
             "boundary_elements",
             *std::max_element(numBndElements, numBndElements + 3),
             &ncDimBndElem);

  // Create netcdf variables
  int ncVarElemSize;
  checkNcError(nc_def_var(ncFile, "element_size", NC_INT, 1, &ncDimPart, &ncVarElemSize));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemSize, NC_COLLECTIVE));

  int ncVarElemVertices;
  int dimsElemVertices[] = {ncDimPart, ncDimElem, ncDimElemVertices};
  checkNcError(
      nc_def_var(ncFile, "element_vertices", NC_INT, 3, dimsElemVertices, &ncVarElemVertices));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemVertices, NC_COLLECTIVE));

  int ncVarElemNeighbors;
  int dimsElemSides[] = {ncDimPart, ncDimElem, ncDimElemSides};
  checkNcError(
      nc_def_var(ncFile, "element_neighbors", NC_INT, 3, dimsElemSides, &ncVarElemNeighbors));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemNeighbors, NC_COLLECTIVE));

  int ncVarElemBoundaries;
  checkNcError(
      nc_def_var(ncFile, "element_boundaries", NC_INT, 3, dimsElemSides, &ncVarElemBoundaries));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemBoundaries, NC_COLLECTIVE));

  int ncVarElemNeighborSides;
  checkNcError(nc_def_var(
      ncFile, "element_neighbor_sides", NC_INT, 3, dimsElemSides, &ncVarElemNeighborSides));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborSides, NC_COLLECTIVE));

  int ncVarElemSideOrientations;
  checkNcError(nc_def_var(
      ncFile, "element_side_orientations", NC_INT, 3, dimsElemSides, &ncVarElemSideOrientations));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemSideOrientations, NC_COLLECTIVE));

  int ncVarElemNeighborRanks;
  checkNcError(nc_def_var(
      ncFile, "element_neighbor_ranks", NC_INT, 3, dimsElemSides, &ncVarElemNeighborRanks));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborRanks, NC_COLLECTIVE));

  int ncVarElemMPIIndices;
  checkNcError(
      nc_def_var(ncFile, "element_mpi_indices", NC_INT, 3, dimsElemSides, &ncVarElemMPIIndices));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemMPIIndices, NC_COLLECTIVE));

  int ncVarElemGroup;
  checkNcError(nc_def_var(ncFile, "element_group", NC_INT, 2, dimsElemSides, &ncVarElemGroup));
  //      checkNcError(nc_var_par_access(ncFile, ncVarElemGroup, NC_COLLECTIVE));

  int ncVarVrtxSize;
  checkNcError(nc_def_var(ncFile, "vertex_size", NC_INT, 1, &ncDimPart, &ncVarVrtxSize));
  //      checkNcError(nc_var_par_access(ncFile, ncVarVrtxSize, NC_COLLECTIVE));

  int ncVarVrtxCoords;
  int dimsVrtxCoords[] = {ncDimPart, ncDimVrtx, ncDimDimension};
  checkNcError(
      nc_def_var(ncFile, "vertex_coordinates", NC_DOUBLE, 3, dimsVrtxCoords, &ncVarVrtxCoords));
  //      checkNcError(nc_var_par_access(ncFile, ncVarVrtxCoords, NC_COLLECTIVE));

  int ncVarBndSize;
  checkNcError(nc_def_var(ncFile, "boundary_size", NC_INT, 1, &ncDimPart, &ncVarBndSize));
  //      checkNcError(nc_var_par_access(ncFile, ncVarBndSize, NC_COLLECTIVE));

  int ncVarBndElemSize;
  int dimsBndElemSize[] = {ncDimPart, ncDimBnd};
  checkNcError(
      nc_def_var(ncFile, "boundary_element_size", NC_INT, 2, dimsBndElemSize, &ncVarBndElemSize));
  //      checkNcError(nc_var_par_access(ncFile, ncVarBndElemSize, NC_COLLECTIVE));

  int ncVarBndElemRank;
  checkNcError(
      nc_def_var(ncFile, "boundary_element_rank", NC_INT, 2, dimsBndElemSize, &ncVarBndElemRank));
  //      checkNcError(nc_var_par_access(ncFile, ncVarBndElemRank, NC_COLLECTIVE));

  int ncVarBndElemLocalIds;
  int dimsBndElemLocalIds[] = {ncDimPart, ncDimBnd, ncDimBndElem};
  checkNcError(nc_def_var(
      ncFile, "boundary_element_localids", NC_INT, 3, dimsBndElemLocalIds, &ncVarBndElemLocalIds));
  //      checkNcError(nc_var_par_access(ncFile, ncVarBndElemLocalIds, NC_COLLECTIVE));

  int writes_done = 0;
  loadBar(writes_done, netcdf_writes);

  // Elements
  int* elemSize = new int[numPartitions[3]];
  std::fill(elemSize, elemSize + numPartitions[3], numElemPerPart[3]);
  checkNcError(nc_put_var_int(ncFile, ncVarElemSize, elemSize));
  delete[] elemSize;
  writes_done++;
  loadBar(writes_done, netcdf_writes);

  std::vector<CubeVertex> vertices;
  vertices.resize(numElemPerPart[3] * 4);
#pragma omp parallel for
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        int odd = (zz + yy + xx) % 2;

        for (unsigned int i = 0; i < 5; i++) {
          for (unsigned int j = 0; j < 4; j++) {
            CubeVertex v;
            v.v[0] = TET_VERTICES[odd][i][j][0] + xx;
            v.v[1] = TET_VERTICES[odd][i][j][1] + yy;
            v.v[2] = TET_VERTICES[odd][i][j][2] + zz;
            vertices[c] = v;
            c++;
          }
        }
      }
    }
  }

  int* elemVertices = new int[numElemPerPart[3] * 4];
  std::map<CubeVertex, int> vertexMap;
  for (unsigned int i = 0; i < vertices.size(); i++) {
    std::map<CubeVertex, int>::iterator it = vertexMap.find(vertices[i]);
    if (it != vertexMap.end()) {
      elemVertices[i] = it->second;
    } else {
      int n = vertexMap.size();
      vertexMap[vertices[i]] = n;
      elemVertices[i] = n;
    }
  }

  for (unsigned int i = 0; i < numPartitions[3]; i++) {
    size_t start[3] = {i, 0, 0};
    size_t count[3] = {1, numElemPerPart[3], 4};
    checkNcError(nc_put_vara_int(ncFile, ncVarElemVertices, start, count, elemVertices));
    writes_done++;
    loadBar(writes_done, netcdf_writes);
  }
  delete[] elemVertices;

  int* elemNeighbors = new int[numElemPerPart[3] * 4];
  const int TET_NEIGHBORS[2][5 * 4] = {
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

#pragma omp parallel for
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        int odd = (zz + yy + xx) % 2;

        memcpy(&elemNeighbors[c], TET_NEIGHBORS[odd], sizeof(int) * 20);
        int offset = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 5;
        for (int i = 0; i < 20; i++)
          elemNeighbors[c + i] += offset;

        if (xx == 0) { // first cube in a partition in x dimension
          if (odd) {
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
          if (odd) {
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
          if (odd) {
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
          if (odd) {
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
          if (odd) {
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
          if (odd) {
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

  for (unsigned int i = 0; i < numPartitions[3]; i++) {
    size_t start[3] = {i, 0, 0};
    size_t count[3] = {1, numElemPerPart[3], 4};
    checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighbors, start, count, elemNeighbors));
    writes_done++;
    loadBar(writes_done, netcdf_writes);
  }
  delete[] elemNeighbors;

  int* elemBoundaries = new int[numElemPerPart[3] * 4];
  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        memset(elemBoundaries, 0, sizeof(int) * numElemPerPart[3] * 4);

        if (x == 0) { // first partition in x dimension
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy) % 2;
              if (odd) {
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

#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy + 1) % 2;
              if (odd) {
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
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx) % 2;
              if (odd) {
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
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx + 1) % 2;
              if (odd) {
                elemBoundaries[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               3] = boundaryMaxy;
                elemBoundaries[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               14] = boundaryMaxy;
              } else {
                elemBoundaries[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               7] = boundaryMaxy;
                elemBoundaries[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               13] = boundaryMaxy;
              }
            }
          }
        }
        if (z == 0) { // first partition in z dimension
#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (yy + xx) % 2;
              if (odd) {
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
#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              //                                                      int odd = (yy+xx+1) % 2;
              //                                                      if (odd) {
              elemBoundaries[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                  numCubesPerPart[0] +
                              xx) *
                                 20 +
                             11] = boundaryMaxz;
              elemBoundaries[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                  numCubesPerPart[0] +
                              xx) *
                                 20 +
                             15] = boundaryMaxz;
              //                                                      } else {
              //                                                              elemBoundaries[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 11] =
              //                                                              boundaryMaxz;
              //                                                              elemBoundaries[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 15] =
              //                                                              boundaryMaxz;
              //                                                      }
            }
          }
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};
        checkNcError(nc_put_vara_int(ncFile, ncVarElemBoundaries, start, count, elemBoundaries));
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }
  delete[] elemBoundaries;

  int* elemNeighborSides = new int[numElemPerPart[3] * 4];
  int* elemNeighborSidesDef = new int[numElemPerPart[3] * 4];
#pragma omp parallel for
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        int odd = (zz + yy + xx) % 2;
        unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        memcpy(&elemNeighborSidesDef[c], TET_SIDE_NEIGHBORS[odd], sizeof(int) * 20);
      }
    }
  }

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        memcpy(elemNeighborSides, elemNeighborSidesDef, sizeof(int) * numElemPerPart[3] * 4);

        if (boundaryMinx != 6 && x == 0) { // first partition in x dimension
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy) % 2;
              if (odd) {
                elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] = 0;
                elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 10] =
                    0;
              } else {
                elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] = 0;
                elemNeighborSides[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 12] =
                    0;
              }
            }
          }
        }
        if (boundaryMaxx != 6 && x == numPartitions[0] - 1) { // last partition in x dimension
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy + 1) % 2;
              if (odd) {
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
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx) % 2;
              if (odd) {
                elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 6] = 0;
                elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 9] = 0;
              } else {
                elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 1] = 0;
                elemNeighborSides[(zz * numCubesPerPart[1] * numCubesPerPart[0] + xx) * 20 + 10] =
                    0;
              }
            }
          }
        }
        if (boundaryMaxy != 6 && y == numPartitions[1] - 1) { // last partition in y dimension
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx + 1) % 2;
              if (odd) {
                elemNeighborSides[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  3] = 0;
                elemNeighborSides[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  14] = 0;
              } else {
                elemNeighborSides[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  7] = 0;
                elemNeighborSides[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  13] = 0;
              }
            }
          }
        }
        if (boundaryMinz != 6 && z == 0) { // first partition in z dimension
#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (yy + xx) % 2;
              if (odd) {
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
#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              //                                                      int odd = (yy+xx+1) % 2;
              //                                                      if (odd) {
              elemNeighborSides[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                     numCubesPerPart[0] +
                                 xx) *
                                    20 +
                                11] = 0;
              elemNeighborSides[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                     numCubesPerPart[0] +
                                 xx) *
                                    20 +
                                15] = 0;
              //                                                      } else {
              //                                                              elemSideNeighbors[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 11] = 0;
              //                                                              elemSideNeighbors[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 15] = 0;
              //                                                      }
            }
          }
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};
        checkNcError(
            nc_put_vara_int(ncFile, ncVarElemNeighborSides, start, count, elemNeighborSides));
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }
  delete[] elemNeighborSides;
  delete[] elemNeighborSidesDef;

  int* elemSideOrientations = new int[numElemPerPart[3] * 4];
  int* elemSideOrientationsDef = new int[numElemPerPart[3] * 4];
#pragma omp parallel for
  for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
    for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
      for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
        int odd = (zz + yy + xx) % 2;
        unsigned int c = ((zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) * 20;
        memcpy(&elemSideOrientationsDef[c], TET_SIDE_ORIENTATIONS[odd], sizeof(int) * 20);
      }
    }
  }

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        memcpy(elemSideOrientations, elemSideOrientationsDef, sizeof(int) * numElemPerPart[3] * 4);

        if (boundaryMinx != 6 && x == 0) { // first partition in x dimension
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy) % 2;
              if (odd) {
                elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20] = 0;
                elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 +
                                     10] = 0;
              } else {
                elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 + 2] =
                    0;
                elemSideOrientations[(zz * numCubesPerPart[1] + yy) * numCubesPerPart[0] * 20 +
                                     12] = 0;
              }
            }
          }
        }
        if (boundaryMaxx != 6 && x == numPartitions[0] - 1) { // last partition in x dimension
#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy + 1) % 2;
              if (odd) {
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
        //                                      #pragma omp parallel for
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
        //                                      #pragma omp parallel for
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
#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (yy + xx) % 2;
              if (odd) {
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
#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              //                                                      int odd = (yy+xx+1) % 2;
              //                                                      if (odd) {
              elemSideOrientations[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                        numCubesPerPart[0] +
                                    xx) *
                                       20 +
                                   11] = 0;
              elemSideOrientations[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                        numCubesPerPart[0] +
                                    xx) *
                                       20 +
                                   15] = 0;
              //                                                      } else {
              //                                                              elemSideOrientations[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 11] = 0;
              //                                                              elemSideOrientations[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 15] = 0;
              //                                                      }
            }
          }
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};
        checkNcError(
            nc_put_vara_int(ncFile, ncVarElemSideOrientations, start, count, elemSideOrientations));
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }
  delete[] elemSideOrientations;
  delete[] elemSideOrientationsDef;

  int* elemNeighborRanks = new int[numElemPerPart[3] * 4];

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        int myrank = (z * numPartitions[1] + y) * numPartitions[0] + x;

        std::fill(elemNeighborRanks, elemNeighborRanks + numElemPerPart[3] * 4, myrank);

        if ((boundaryMinx == 6 && numPartitions[0] > 1) ||
            x != 0) { // first partition in x dimension
          int rank = (z * numPartitions[1] + y) * numPartitions[0] +
                     (x - 1 + numPartitions[0]) % numPartitions[0];

#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy) % 2;
              if (odd) {
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
          int rank = (z * numPartitions[1] + y) * numPartitions[0] + (x + 1) % numPartitions[0];

#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy + 1) % 2;
              if (odd) {
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
        if ((boundaryMiny == 6 && numPartitions[1] > 1) ||
            y != 0) { // first partition in y dimension
          int rank = (z * numPartitions[1] + (y - 1 + numPartitions[1]) % numPartitions[1]) *
                         numPartitions[0] +
                     x;

#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx) % 2;
              if (odd) {
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
          int rank = (z * numPartitions[1] + (y + 1) % numPartitions[1]) * numPartitions[0] + x;

#pragma omp parallel for
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx + 1) % 2;
              if (odd) {
                elemNeighborRanks[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  3] = rank;
                elemNeighborRanks[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  14] = rank;
              } else {
                elemNeighborRanks[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  7] = rank;
                elemNeighborRanks[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                       numCubesPerPart[0] +
                                   xx) *
                                      20 +
                                  13] = rank;
              }
            }
          }
        }
        if ((boundaryMinz == 6 && numPartitions[2] > 1) ||
            z != 0) { // first partition in z dimension
          int rank = (((z - 1 + numPartitions[2]) % numPartitions[2]) * numPartitions[1] + y) *
                         numPartitions[0] +
                     x;

#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (yy + xx) % 2;
              if (odd) {
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
          int rank = (((z + 1) % numPartitions[2]) * numPartitions[1] + y) * numPartitions[0] + x;

#pragma omp parallel for
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              //                                                      int odd = (yy+xx+1) % 2;
              //                                                      if (odd) {
              elemNeighborRanks[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                     numCubesPerPart[0] +
                                 xx) *
                                    20 +
                                11] = rank;
              elemNeighborRanks[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                     numCubesPerPart[0] +
                                 xx) *
                                    20 +
                                15] = rank;
              //                                                      } else {
              //                                                              elemNeighborRanks[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 11] = rank;
              //                                                              elemNeighborRanks[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
              //                                                              * 20 + 15] = rank;
              //                                                      }
            }
          }
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};
        checkNcError(
            nc_put_vara_int(ncFile, ncVarElemNeighborRanks, start, count, elemNeighborRanks));
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }
  delete[] elemNeighborRanks;

  int* elemMPIIndices = new int[numElemPerPart[3] * 4];
  int* bndLocalIds = new int[*std::max_element(numBndElements, numBndElements + 3)];

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        memset(elemMPIIndices, 0, sizeof(int) * numElemPerPart[3] * 4);

        unsigned int bndSize = 0;

        if ((boundaryMinz == 6 && numPartitions[2] > 1) ||
            z != 0) { // first partition in z dimension
          int nextMPIIndex = 0;
          for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (yy + xx) % 2;
              if (odd) {
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

          size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0u};
          size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          int rank = (((z - 1 + numPartitions[2]) % numPartitions[2]) * numPartitions[1] + y) *
                         numPartitions[0] +
                     x;
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));

          bndSize++;
        }
        if ((boundaryMiny == 6 && numPartitions[1] > 1) ||
            y != 0) { // first partition in y dimension
          int nextMPIIndex = 0;
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx) % 2;
              if (odd) {
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

          size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0u};
          size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          int rank = (z * numPartitions[1] + (y - 1 + numPartitions[1]) % numPartitions[1]) *
                         numPartitions[0] +
                     x;
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));

          bndSize++;
        }
        if ((boundaryMinx == 6 && numPartitions[0] > 1) ||
            x != 0) { // first partition in x dimension
          int nextMPIIndex = 0;
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy) % 2;
              if (odd) {
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

          size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0u};
          size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          int rank = (z * numPartitions[1] + y) * numPartitions[0] +
                     (x - 1 + numPartitions[0]) % numPartitions[0];
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));

          bndSize++;
        }
        if ((boundaryMaxx == 6 && numPartitions[0] > 1) ||
            x != numPartitions[0] - 1) { // last partition in x dimension
          int nextMPIIndex = 0;
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int yy = 0; yy < numCubesPerPart[1]; yy++) {
              int odd = (zz + yy + 1) % 2;
              if (odd) {
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

          size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0};
          size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          int rank = (z * numPartitions[1] + y) * numPartitions[0] + (x + 1) % numPartitions[0];
          rank = (rank + numPartitions[3]) % numPartitions[3];
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));

          bndSize++;
        }
        if ((boundaryMaxy == 6 && numPartitions[1] > 1) ||
            y != numPartitions[1] - 1) { // last partition in y dimension
          int nextMPIIndex = 0;
          for (unsigned int zz = 0; zz < numCubesPerPart[2]; zz++) {
            for (unsigned int xx = 0; xx < numCubesPerPart[0]; xx++) {
              int odd = (zz + xx + 1) % 2;
              if (odd) {
                bndLocalIds[nextMPIIndex] =
                    ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                    5;
                elemMPIIndices[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               3] = nextMPIIndex++;
                bndLocalIds[nextMPIIndex] =
                    ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                        5 +
                    3;
                elemMPIIndices[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               14] = nextMPIIndex++;
              } else {
                bndLocalIds[nextMPIIndex] =
                    ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                        5 +
                    1;
                elemMPIIndices[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               7] = nextMPIIndex++;
                bndLocalIds[nextMPIIndex] =
                    ((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) * numCubesPerPart[0] + xx) *
                        5 +
                    3;
                elemMPIIndices[((zz * numCubesPerPart[1] + numCubesPerPart[1] - 1) *
                                    numCubesPerPart[0] +
                                xx) *
                                   20 +
                               13] = nextMPIIndex++;
              }
            }
          }

          size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0};
          size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          int rank = (z * numPartitions[1] + (y + 1) % numPartitions[1]) * numPartitions[0] + x;
          rank = (rank + numPartitions[3]) % numPartitions[3];
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));

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
              elemMPIIndices[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                  numCubesPerPart[0] +
                              xx) *
                                 20 +
                             11] = nextMPIIndex++;
              bndLocalIds[nextMPIIndex] =
                  (((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) * numCubesPerPart[0] + xx) *
                      5 +
                  3;
              elemMPIIndices[(((numCubesPerPart[2] - 1) * numCubesPerPart[1] + yy) *
                                  numCubesPerPart[0] +
                              xx) *
                                 20 +
                             15] = nextMPIIndex++;
            }
          }

          size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, bndSize, 0u};
          size_t count[3] = {1, 1, static_cast<unsigned int>(nextMPIIndex)};
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          int rank = (((z + 1) % numPartitions[2]) * numPartitions[1] + y) * numPartitions[0] + x;
          rank = (rank + numPartitions[3]) % numPartitions[3];
          checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));

          bndSize++;
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};
        checkNcError(nc_put_vara_int(ncFile, ncVarElemMPIIndices, start, count, elemMPIIndices));

        checkNcError(nc_put_var1_uint(ncFile, ncVarBndSize, &start[0], &bndSize));
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }
  delete[] elemMPIIndices;
  delete[] bndLocalIds;

  // Set material zone to 1
  int* elemGroup = new int[numElemPerPart[3]];
  std::fill(elemGroup, elemGroup + numElemPerPart[3], 1);
  for (unsigned int x = 0; x < numPartitions[3]; x++) {
    size_t start[2] = {x, 0};
    size_t count[2] = {1, numElemPerPart[3]};
    checkNcError(nc_put_vara_int(ncFile, ncVarElemGroup, start, count, elemGroup));
  }
  delete[] elemGroup;

  // Vertices
  std::map<int, CubeVertex> uniqueVertices;
  transform(vertexMap.begin(),
            vertexMap.end(),
            inserter(uniqueVertices, uniqueVertices.begin()),
            flip_pair<CubeVertex, int>);

  int* vrtxSize = new int[numPartitions[3]];
  std::fill(vrtxSize, vrtxSize + numPartitions[3], uniqueVertices.size());
  checkNcError(nc_put_var_int(ncFile, ncVarVrtxSize, vrtxSize));
  delete[] vrtxSize;
  writes_done++;
  loadBar(writes_done, netcdf_writes);

  double* vrtxCoords = new double[uniqueVertices.size() * 3];

  double halfWidthX = scaleX / 2.0;
  double halfWidthY = scaleY / 2.0;
  double halfWidthZ = scaleZ / 2.0;

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      for (unsigned int x = 0; x < numPartitions[0]; x++) {

#pragma omp parallel for
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

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, uniqueVertices.size(), 3};
        checkNcError(nc_put_vara_double(ncFile, ncVarVrtxCoords, start, count, vrtxCoords));
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }

  delete[] vrtxCoords;

  checkNcError(nc_close(ncFile));

  logInfo() << "Finished";
}
#endif // USE_NETCDF

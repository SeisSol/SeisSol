#include "utils/logger.h"
#include "utils/env.h"
#include "utils/args.h"
#include "CubeGenerator.h"
#include "MeshReader.h"

#ifdef USE_NETCDF
#ifdef USE_MPI
#include "Parallel/MPI.h"
#endif // USE_MPI

#include <netcdf.h>
#include <netcdf_par.h>

#include <omp.h>

#include <algorithm>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

#include <iostream>
#include <sstream>
#include <string>

namespace {
using t_vertex = int[3];

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
} // anonymous namespace




seissol::geometry::CubeGenerator::CubeGenerator(
    int rank,
    int nProcs,
    const std::string& meshFile,
    const seissol::geometry::CubeGeneratorParameters& cubeParams)
    : seissol::geometry::MeshReader(rank), // init base class
      rank(rank), nProcs(nProcs) {
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
    // check if numCubes is multiple of numPartitions, should only fail in numPartitions[0]
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
                                                     const std::string& meshFile) {

  logInfo() << "Total number of cubes:" << numCubes[0] << 'x' << numCubes[1] << 'x' << numCubes[2]
            << '=' << numCubes[3];
  logInfo() << "Total number of partitions" << numPartitions[0] << 'x' << numPartitions[1] << 'x'
            << numPartitions[2] << '=' << numPartitions[3];
  logInfo() << "Total number of cubes per partition:" << numCubesPerPart[0] << 'x'
            << numCubesPerPart[1] << 'x' << numCubesPerPart[2] << '=' << numCubesPerPart[3];
  //TODO: do 5 * numCubesPerPart[0-2]
  logInfo() << "Total number of elements per partition:" << numElemPerPart[0] << 'x'
            << numElemPerPart[1] << 'x' << numElemPerPart[2] << '=' << numElemPerPart[3];
  logInfo() << "Using" << omp_get_max_threads() << "threads";

  int netcdf_writes = 2 + numPartitions[3] * 8;

  // Create the netcdf file
  int ncFile;
  int masterRank;
  unsigned int groupSize = 1;
  // TODO: MPI is technically not needed for the CubeGenerator, groupSize not needed as well
#ifdef USE_MPI
  groupSize = utils::Env::get<unsigned int>("SEISSOL_NETCDF_GROUP_SIZE", 1);
  if (nProcs % groupSize != 0)
    logError() << "#Processes must be a multiple of the group size" << groupSize;

  int master = (rank / groupSize) * groupSize;
  MPI_Comm commMaster;
  MPI_Comm_split(
      seissol::MPI::mpi.comm(), rank % groupSize == 0 ? 1 : MPI_UNDEFINED, rank, &commMaster);

  masterRank = -1;

  if (commMaster != MPI_COMM_NULL) {
    MPI_Comm_rank(commMaster, &masterRank);
    checkNcError(nc_open_par(meshFile.c_str(), NC_NETCDF4 | NC_MPIIO, commMaster, MPI_INFO_NULL, &ncFile));
  }
#else  // USE_MPI
  masterRank = rank; // = 0;
  checkNcError(nc_open(meshFile, NC_NETCDF4, &ncFile));
#endif // USE_MPI
       /* TODO: these variables are used to define properties in the .nc file
                redundant as we now read a finished .nc file instead of creating one
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
       */

  size_t bndSize = -1; // 
  size_t bndElemSize = -1;
  
//  int ncVarElemSize = -1; //TODO: numElemPerPart[3] * numPartitions[3]
//  int ncVarElemVertices = -1;
//  int ncVarElemNeighbors = -1;
//  int ncVarElemBoundaries = -1;
//  int ncVarElemNeighborSides = -1;
//  int ncVarElemSideOrientations = -1;
//  int ncVarElemNeighborRanks = -1;
//  int ncVarElemMPIIndices = -1;
//  int ncVarElemGroup = -1;
  bool hasGroup = false;
//  int ncVarVrtxSize = -1;
//  int ncVarVrtxCoords = -1;
//  int ncVarBndSize = -1; // TODO: bndSize gets calculated anyway
  int ncVarBndElemSize = -1;
  int ncVarBndElemRank = -1;
  int ncVarBndElemLocalIds = -1;

  int* sizes = 0L;
  int maxSize = 0;
  if (masterRank == 0) {
    // Get important dimensions
    //int ncDimPart;
    //checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
    size_t partitions = numPartitions[3]; // TODO: dimlen = numPartitions[3] get numPartitions
    //checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &partitions));

    if (partitions != static_cast<unsigned int>(nProcs))
      logError() << "Number of partitions in netCDF file does not match number of MPI ranks.";

    //int ncDimBndSize; // TODO: acc. to CubeGenerator: "boundaries" is 6 -> bndSize is also 6: from netcdf definiton:
    // 3rd parameter in nc_def_dim is the length of the dimension to be created!
    // TODO: nic_inq_dimlen returns the const value that the CubeGenerator would define using nc_def_dim!
    /*checkNcError(nc_inq_dimid(ncFile, "boundaries", &ncDimBndSize));
    checkNcError(nc_inq_dimlen(ncFile, ncDimBndSize, &bndSize));*/
    bndSize = 6;

    //int ncDimBndElem; // TODO: bndElemSize nc_inq_dimlen = *std::max_element(numBndElements, numBndElements+3)
    /*checkNcError(nc_inq_dimid(ncFile, "boundary_elements", &ncDimBndElem));
    checkNcError(nc_inq_dimlen(ncFile, ncDimBndElem, &bndElemSize));*/
    bndElemSize = *std::max_element(numBndElements, numBndElements+3);
  }

#ifdef USE_MPI
  // Broadcast boundary sizes
  unsigned long buf[2] = {bndSize, bndElemSize};
  MPI_Bcast(buf, 2, MPI_UNSIGNED_LONG, 0, seissol::MPI::mpi.comm());
  bndSize = buf[0];
  bndElemSize = buf[1];
#endif // USE_MPI

  // Create netcdf variables
  if (masterRank >= 0) {
/*    //    int ncVarElemSize;
    //      checkNcError(nc_def_var(ncFile, "element_size", NC_INT, 1, &ncDimPart, &ncVarElemSize));
    checkNcError(nc_inq_varid(ncFile, "element_size", &ncVarElemSize));
    collectiveAccess(ncFile, ncVarElemSize);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemSize, NC_COLLECTIVE));

    //      int dimsElemVertices[] = {ncDimPart, ncDimElem, ncDimElemVertices};
    //    int ncVarElemVertices;
    //      checkNcError(
    //          nc_def_var(ncFile, "element_vertices", NC_INT, 3, dimsElemVertices,
    //          &ncVarElemVertices));
    checkNcError(nc_inq_varid(ncFile, "element_vertices", &ncVarElemVertices));
    collectiveAccess(ncFile, ncVarElemVertices);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemVertices, NC_COLLECTIVE));

    //      int dimsElemSides[] = {ncDimPart, ncDimElem, ncDimElemSides};
    //    int ncVarElemNeighbors;
    //      checkNcError(
    //          nc_def_var(ncFile, "element_neighbors", NC_INT, 3, dimsElemSides,
    //          &ncVarElemNeighbors));
    checkNcError(nc_inq_varid(ncFile, "element_neighbors", &ncVarElemNeighbors));
    collectiveAccess(ncFile, ncVarElemNeighbors);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemNeighbors, NC_COLLECTIVE));

    //    int ncVarElemBoundaries;
    //      checkNcError(
    //          nc_def_var(ncFile, "element_boundaries", NC_INT, 3, dimsElemSides,
    //          &ncVarElemBoundaries));
    checkNcError(nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries));
    collectiveAccess(ncFile, ncVarElemBoundaries);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemBoundaries, NC_COLLECTIVE));

    //    int ncVarElemNeighborSides;
    //      checkNcError(nc_def_var(
    //          ncFile, "element_neighbor_sides", NC_INT, 3, dimsElemSides,
    //          &ncVarElemNeighborSides));
    checkNcError(nc_inq_varid(ncFile, "element_neighbor_sides", &ncVarElemNeighborSides));
    collectiveAccess(ncFile, ncVarElemNeighborSides);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborSides, NC_COLLECTIVE));

    //    int ncVarElemSideOrientations;
    //      checkNcError(nc_def_var(
    //          ncFile, "element_side_orientations", NC_INT, 3, dimsElemSides,
    //          &ncVarElemSideOrientations));
    checkNcError(nc_inq_varid(ncFile, "element_side_orientations", &ncVarElemSideOrientations));
    collectiveAccess(ncFile, ncVarElemSideOrientations);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemSideOrientations, NC_COLLECTIVE));

    //    int ncVarElemNeighborRanks;
    //      checkNcError(nc_def_var(
    //          ncFile, "element_neighbor_ranks", NC_INT, 3, dimsElemSides,
    //          &ncVarElemNeighborRanks));
    checkNcError(nc_inq_varid(ncFile, "element_neighbor_ranks", &ncVarElemNeighborRanks));
    collectiveAccess(ncFile, ncVarElemNeighborRanks);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborRanks, NC_COLLECTIVE));

    //    int ncVarElemMPIIndices;
    //      checkNcError(
    //          nc_def_var(ncFile, "element_mpi_indices", NC_INT, 3, dimsElemSides,
    //          &ncVarElemMPIIndices));
    checkNcError(nc_inq_varid(ncFile, "element_mpi_indices", &ncVarElemMPIIndices));
    collectiveAccess(ncFile, ncVarElemMPIIndices);
    //      checkNcError(nc_var_par_access(ncFile, ncVarElemMPIIndices, NC_COLLECTIVE));

    //    int ncVarElemGroup;
    //      checkNcError(nc_def_var(ncFile, "element_group", NC_INT, 2, dimsElemSides,
    //      &ncVarElemGroup));
    //      //      checkNcError(nc_var_par_access(ncFile, ncVarElemGroup, NC_COLLECTIVE));
    int ncResult = nc_inq_varid(ncFile, "element_group", &ncVarElemGroup);
    if (ncResult != NC_ENOTVAR) {  // TODO: why is this the only variable handled like this?
      checkNcError(ncResult);
      hasGroup = true; // TODO: necessary for CubeGenerator?
      collectiveAccess(ncFile, ncVarElemGroup);
    }

    //    int ncVarVrtxSize;
    //      checkNcError(nc_def_var(ncFile, "vertex_size", NC_INT, 1, &ncDimPart, &ncVarVrtxSize));
    checkNcError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
    collectiveAccess(ncFile, ncVarVrtxSize);
    //      checkNcError(nc_var_par_access(ncFile, ncVarVrtxSize, NC_COLLECTIVE));

    //      int dimsVrtxCoords[] = {ncDimPart, ncDimVrtx, ncDimDimension};
    //    int ncVarVrtxCoords;
    //      checkNcError(
    //           nc_def_var(ncFile, "vertex_coordinates", NC_DOUBLE, 3, dimsVrtxCoords,
    //           &ncVarVrtxCoords));
    checkNcError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
    collectiveAccess(ncFile, ncVarVrtxCoords);
    //      checkNcError(nc_var_par_access(ncFile, ncVarVrtxCoords, NC_COLLECTIVE));
*/
    //    int ncVarBndSize;
    //      checkNcError(nc_def_var(ncFile, "boundary_size", NC_INT, 1, &ncDimPart, &ncVarBndSize));
    //checkNcError(nc_inq_varid(ncFile, "boundary_size", &ncVarBndSize));
    //collectiveAccess(ncFile, ncVarBndSize);
    //      checkNcError(nc_var_par_access(ncFile, ncVarBndSize, NC_COLLECTIVE));

    //      int dimsBndElemSize[] = {ncDimPart, ncDimBnd};
    //    int ncVarBndElemSize;
    //      checkNcError(
    //          nc_def_var(ncFile, "boundary_element_size", NC_INT, 2, dimsBndElemSize,
    //          &ncVarBndElemSize));
    checkNcError(nc_inq_varid(ncFile, "boundary_element_size", &ncVarBndElemSize));
    collectiveAccess(ncFile, ncVarBndElemSize);
    //      checkNcError(nc_var_par_access(ncFile, ncVarBndElemSize, NC_COLLECTIVE));

    //    int ncVarBndElemRank;
    //      checkNcError(
    //          nc_def_var(ncFile, "boundary_element_rank", NC_INT, 2, dimsBndElemSize,
    //          &ncVarBndElemRank));
    checkNcError(nc_inq_varid(ncFile, "boundary_element_rank", &ncVarBndElemRank));
    collectiveAccess(ncFile, ncVarBndElemRank);
    //      checkNcError(nc_var_par_access(ncFile, ncVarBndElemRank, NC_COLLECTIVE));

    //      int dimsBndElemLocalIds[] = {ncDimPart, ncDimBnd, ncDimBndElem};
    //    int ncVarBndElemLocalIds;
    //      checkNcError(nc_def_var(
    //          ncFile, "boundary_element_localids", NC_INT, 3, dimsBndElemLocalIds,
    //          &ncVarBndElemLocalIds));
    checkNcError(nc_inq_varid(ncFile, "boundary_element_localids", &ncVarBndElemLocalIds));
    collectiveAccess(ncFile, ncVarBndElemLocalIds);
    //      checkNcError(nc_var_par_access(ncFile, ncVarBndElemLocalIds, NC_COLLECTIVE));

    // TODO: reading of variables is done, below we start with the elements
    // Elements
    sizes = new int[groupSize];

    for (int i = groupSize - 1; i >= 0; i--) {
      size_t start = static_cast<size_t>(i + rank);

      // TODO: remove numPartitions[3] -> every MPI process calls CubeGenerator with its corresponding rank
      int size = numElemPerPart[3]; // * numPartitions[3]; // TODO: see nc var decls
      //checkNcError(nc_get_var1_int(ncFile, ncVarElemSize, &start, &size));

      sizes[i] = size;
      if (i != 0) {
#ifdef USE_MPI
        // Forward size (except the own)
        MPI_Send(&size, 1, MPI_INT, i + rank, 0, seissol::MPI::mpi.comm());
#else  // USE_MPI
        assert(false);
#endif // USE_MPI
      }

      maxSize = std::max(maxSize, size);
    }
// TODO: delete #endif
  } else {
#ifdef USE_MPI
    sizes = new int[1];
    MPI_Recv(sizes, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);

    maxSize = sizes[0];
#else  // USE_MPI
  assert(false);
#endif // USE_MPI
  }

#ifdef USE_MPI
  // Broadcast group information
  // TODO: delete?
  int iHasGroup = hasGroup;
  MPI_Bcast(&iHasGroup, 1, MPI_INT, 0, seissol::MPI::mpi.comm());
  hasGroup = iHasGroup != 0;
#endif // USE_MPI

  int writes_done = 0;
  loadBar(writes_done, netcdf_writes);

  m_elements.resize(sizes[0]);

  // TODO: delete?
/*  int* elemSize = new int[numPartitions[3]];
  std::fill(elemSize, elemSize + numPartitions[3], numElemPerPart[3]);
  // checkNcError(nc_put_var_int(ncFile, ncVarElemSize, elemSize));
  delete[] elemSize;*/
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

/*  for (unsigned int i = 0; i < numPartitions[3]; i++) {
    size_t start[3] = {i, 0, 0};
    size_t count[3] = {1, numElemPerPart[3], 4};
    //checkNcError(nc_put_vara_int(ncFile, ncVarElemVertices, start, count, elemVertices));
    writes_done++;
    loadBar(writes_done, netcdf_writes);
  }*/

  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);

//  delete[] elemVertices;

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
    //checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighbors, start, count, elemNeighbors));
    //writes_done++;
    //loadBar(writes_done, netcdf_writes);
  }

  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);
//  delete[] elemNeighbors;

  int* elemBoundaries = new int[numElemPerPart[3] * 4];
  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
      // TODO: change this back if it doesnt work!
//      for (unsigned int x = rank; x <= rank; x++) {
//      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        unsigned int x = rank;
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
                                                                    int odd = (yy+xx+1) % 2;
                                                                    if (odd) {
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
                                                                    } else {
                                                                            elemBoundaries[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
                                                                            * 20 + 11] =
                                                                            boundaryMaxz;
                                                                            elemBoundaries[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
                                                                            * 20 + 15] =
                                                                            boundaryMaxz;
                                                                    }
            }
          }
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};

        ElemNeighbors* elemBoundariesCast = reinterpret_cast<ElemBoundaries*>(elemBoundaries); 
        for (int i = 0; i < sizes[0]; i++) {
          memcpy(m_elements[i].boundaries, &elemBoundariesCast[i], sizeof(ElemBoundaries));
        }

        //checkNcError(nc_put_vara_int(ncFile, ncVarElemBoundaries, start, count, elemBoundaries));
        //writes_done++;
        //loadBar(writes_done, netcdf_writes);
//      }
    }
  }
  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);

//  delete[] elemBoundaries;

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
        //checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborSides, start, count, elemNeighborSides));
        //writes_done++;
        //loadBar(writes_done, netcdf_writes);
      }
    }
  }
  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);

//  delete[] elemNeighborSides;
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
        //checkNcError(nc_put_vara_int(ncFile, ncVarElemSideOrientations, start, count, elemSideOrientations));
        //writes_done++;
        //loadBar(writes_done, netcdf_writes);
      }
    }
  }

  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);

//  delete[] elemSideOrientations;
  delete[] elemSideOrientationsDef;

  int* elemNeighborRanks = new int[numElemPerPart[3] * 4];

  for (unsigned int z = 0; z < numPartitions[2]; z++) {
    for (unsigned int y = 0; y < numPartitions[1]; y++) {
//TODO: change this back if it doesnt work!
//      for (unsigned int x = rank; x <= rank; x++) {
//      for (unsigned int x = 0; x < numPartitions[0]; x++) {
        unsigned int x = rank;
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
//                                                                    int odd = (yy+xx+1) % 2;
//                                                                    if (odd) {
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
//                                                                    } else {
//                                                                            elemNeighborRanks[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
//                                                                            * 20 + 11] = rank;
//                                                                            elemNeighborRanks[(((numCubesPerPart[2]-1)*numCubesPerPart[1]+yy)*numCubesPerPart[0]+xx)
//                                                                            * 20 + 15] = rank;
//                                                                    }
            }
          }
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};

        ElemNeighborRanks* elemNeighborRanksCast = reinterpret_cast<ElemNeighborRanks*>(elemNeighborRanks); 
        for (int i = 0; i < sizes[0]; i++) {
          memcpy(m_elements[i].neighborRanks, &elemNeighborRanksCast[i], sizeof(ElemNeighborRanks));
        }

        //checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborRanks, start, count, elemNeighborRanks));
        //writes_done++;
        //loadBar(writes_done, netcdf_writes);
//      }
    }
  }
  
  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);

//  delete[] elemNeighborRanks;

  int* elemMPIIndices = new int[numElemPerPart[3] * 4];
  int* bndLocalIds = new int[*std::max_element(numBndElements, numBndElements + 3)];

  // TODO: set bndSize outside of loop to 0 so we can access the result at a later stage
  size_t* bndSizePtr = new size_t[numPartitions[3]];

  //TODO: int or size_t?
  int* bndElemSizePtr = new int[numPartitions[3] * bndSize];
  int* bndElemRankPtr = new int[numPartitions[3] * bndSize];
  int* bndElemLocalIdsPtr = new int[numPartitions[3] * bndSize * bndElemSize];
  int* elemMPIIndicesPtr = new int[numPartitions[3] * numElemPerPart[3] * 4];

  int bndSizeGlobal = bndSize;

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
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          //bndElemSizePtr[start[0]][start[1]] = nextMPIIndex;
          int rank = (((z - 1 + numPartitions[2]) % numPartitions[2]) * numPartitions[1] + y) *
                         numPartitions[0] +
                     x;
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          // checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));
          //bndElemRankPtr[start[0]][start[1]] = rank;

          // bndElemSizePtr[numPartitions[3] * bndSize]
          bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
          // bndElemRankPtr[numPartitions[3] * bndSize]
          bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
          // bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
          memcpy(&bndElemLocalIdsPtr[(start[0]*bndSizeGlobal + start[1])*bndElemSize + start[2]], bndLocalIds,
                 sizeof(int)*count[2]);

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
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          //bndElemSizePtr[start[0]][start[1]] = nextMPIIndex;
          int rank = (z * numPartitions[1] + (y - 1 + numPartitions[1]) % numPartitions[1]) *
                         numPartitions[0] +
                     x;
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          // checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));
          //bndElemRankPtr[start[0]][start[1]] = rank;
          //memcpy(&bndElemLocalIdsPtr[start[0]][start[1]], bndLocalIds, sizeof(bndLocalIds) * count[2]);

          // bndElemSizePtr[numPartitions[3] * bndSize]
          bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
          // bndElemRankPtr[numPartitions[3] * bndSize]
          bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
          // bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
          memcpy(&bndElemLocalIdsPtr[(start[0]*bndSizeGlobal + start[1])*bndElemSize + start[2]], bndLocalIds,
                 sizeof(int)*count[2]);

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
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex)); 
          //bndElemSizePtr[start[0]][start[1]] = nextMPIIndex;
          int rank = (z * numPartitions[1] + y) * numPartitions[0] +
                     (x - 1 + numPartitions[0]) % numPartitions[0];
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          // checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));
          //bndElemRankPtr[start[0]][start[1]] = rank;
          //memcpy(&bndElemLocalIdsPtr[start[0]][start[1]], bndLocalIds, sizeof(bndLocalIds) * count[2]);

          // bndElemSizePtr[numPartitions[3] * bndSize]
          bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
          // bndElemRankPtr[numPartitions[3] * bndSize]
          bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
          // bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
          memcpy(&bndElemLocalIdsPtr[(start[0]*bndSizeGlobal + start[1])*bndElemSize + start[2]], bndLocalIds,
                 sizeof(int)*count[2]);

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
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          //bndElemSizePtr[start[0]][start[1]] = nextMPIIndex;
          int rank = (z * numPartitions[1] + y) * numPartitions[0] + (x + 1) % numPartitions[0];
          rank = (rank + numPartitions[3]) % numPartitions[3];
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          // checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));
          //bndElemRankPtr[start[0]][start[1]] = rank;
          //memcpy(&bndElemLocalIdsPtr[start[0]][start[1]], bndLocalIds, sizeof(bndLocalIds) * count[2]);

          // bndElemSizePtr[numPartitions[3] * bndSize]
          bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
          // bndElemRankPtr[numPartitions[3] * bndSize]
          bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
          // bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
          memcpy(&bndElemLocalIdsPtr[(start[0]*bndSizeGlobal + start[1])*bndElemSize + start[2]], bndLocalIds,
                 sizeof(int)*count[2]);

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
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          //bndElemSizePtr[start[0]][start[1]] = nextMPIIndex;
          int rank = (z * numPartitions[1] + (y + 1) % numPartitions[1]) * numPartitions[0] + x;
          rank = (rank + numPartitions[3]) % numPartitions[3];
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          // checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));
          //bndElemRankPtr[start[0]][start[1]] = rank;
          //memcpy(&bndElemLocalIdsPtr[start[0]][start[1]], bndLocalIds, sizeof(bndLocalIds) * count[2]);

          // bndElemSizePtr[numPartitions[3] * bndSize]
          bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
          // bndElemRankPtr[numPartitions[3] * bndSize]
          bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
          // bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
          memcpy(&bndElemLocalIdsPtr[(start[0]*bndSizeGlobal + start[1])*bndElemSize + start[2]], bndLocalIds,
                 sizeof(int)*count[2]);

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
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, start, &nextMPIIndex));
          //bndElemSizePtr[start[0]][start[1]] = nextMPIIndex;
          int rank = (((z + 1) % numPartitions[2]) * numPartitions[1] + y) * numPartitions[0] + x;
          rank = (rank + numPartitions[3]) % numPartitions[3];
          // checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, start, &rank));
          // checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, bndLocalIds));
          //bndElemRankPtr[start[0]][start[1]] = rank;
          //memcpy(&bndElemLocalIdsPtr[start[0]][start[1]], bndLocalIds, sizeof(bndLocalIds) * count[2]);

          // bndElemSizePtr[numPartitions[3] * bndSize]
          bndElemSizePtr[start[0] * bndSizeGlobal + start[1]] = nextMPIIndex;
          // bndElemRankPtr[numPartitions[3] * bndSize]
          bndElemRankPtr[start[0] * bndSizeGlobal + start[1]] = rank;
          // bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
          memcpy(&bndElemLocalIdsPtr[(start[0]*bndSizeGlobal + start[1])*bndElemSize + start[2]], bndLocalIds,
                 sizeof(int)*count[2]);

          bndSize++;
        }

        size_t start[3] = {(z * numPartitions[1] + y) * numPartitions[0] + x, 0, 0};
        size_t count[3] = {1, numElemPerPart[3], 4};
        // checkNcError(nc_put_vara_int(ncFile, ncVarElemMPIIndices, start, count, elemMPIIndices));
        //for (int k = 0; k < count[1]; k++) {
        //  memcpy(elemMPIIndicesPtr[start[0] + k], elemMPIIndices, sizeof(elemMPIIndices) * count[2]);
        //}

        //elemMPIIndicesPtr[numElemPerPart[3] * 4] ->  numElemPerPart[i] = numCubesPerPart[i] * 5 -> numCubes[i] / numPartitions[i] * 5
        //memcpy(&elemMPIIndicesPtr[start[0]], elemMPIIndices, sizeof(*elemMPIIndices)/**count[1]*count[2]*/);
        
        // checkNcError(nc_put_var1_uint(ncFile, ncVarBndSize, &start[0], &bndSize));
        bndSizePtr[start[0]] = bndSize; 
        //writes_done++;
        //loadBar(writes_done, netcdf_writes);
      }
    }
  }

  //TODO: elemMPIIndices    = new int[                 numElemPerPart[3]*4]
  //      elemMPIIndicesPtr = new int[numPartitions[3]*numElemPerPart[3]*4]

  writes_done += numPartitions[3];
  loadBar(writes_done, netcdf_writes);
  //delete[] elemMPIIndices;
  //delete[] bndLocalIds;

  // Set material zone to 1
  int* elemGroup = new int[numElemPerPart[3]];
  std::fill(elemGroup, elemGroup + numElemPerPart[3], 1);
  for (unsigned int x = 0; x < numPartitions[3]; x++) {
    size_t start[2] = {x, 0};
    size_t count[2] = {1, numElemPerPart[3]};
    // checkNcError(nc_put_vara_int(ncFile, ncVarElemGroup, start, count, elemGroup));
  }
//  delete[] elemGroup;

  // TODO: use smart pointers to avoid delete[] calls? if this doesnt work, switch back to pointers
  ElemVertices* elemVerticesCast = reinterpret_cast<ElemVertices*>(elemVertices);
  ElemNeighbors* elemNeighborsCast = reinterpret_cast<ElemNeighbors*>(elemNeighbors); 
  ElemNeighborSides* elemNeighborSidesCast = reinterpret_cast<ElemNeighborSides*>(elemNeighborSides); 
  ElemSideOrientations* elemSideOrientationsCast = reinterpret_cast<ElemSideOrientations*>(elemSideOrientations); 
  ElemBoundaries* elemBoundariesCast = reinterpret_cast<ElemBoundaries*>(elemBoundaries); 
  ElemNeighborRanks* elemNeighborRanksCast = reinterpret_cast<ElemNeighborRanks*>(elemNeighborRanks); 
  ElemMPIIndices* elemMPIIndicesCast = reinterpret_cast<ElemMPIIndices*>(elemMPIIndices); 
  ElemGroup* elemGroupCast = reinterpret_cast<ElemGroup*>(elemGroup); 

  for (int i = 0; i < sizes[0]; i++) {
    m_elements[i].localId = i;

    memcpy(m_elements[i].vertices, &elemVerticesCast[i], sizeof(ElemVertices));
    memcpy(m_elements[i].neighbors, &elemNeighborsCast[i], sizeof(ElemNeighbors));
    memcpy(m_elements[i].neighborSides, &elemNeighborSidesCast[i], sizeof(ElemNeighborSides));
    memcpy(m_elements[i].sideOrientations, &elemSideOrientationsCast[i], sizeof(ElemSideOrientations));
    // TODO: the following 2 memcpys are done in their respective loops above
//    memcpy(m_elements[i].boundaries, &elemBoundariesCast[i], sizeof(ElemBoundaries));
//    memcpy(m_elements[i].neighborRanks, &elemNeighborRanksCast[i], sizeof(ElemNeighborRanks));
    memcpy(m_elements[i].mpiIndices, &elemMPIIndicesCast[i], sizeof(ElemMPIIndices));
    m_elements[i].group = elemGroupCast[i];
 }

  delete[] elemVerticesCast;
  delete[] elemNeighborsCast;
  delete[] elemNeighborSidesCast;
  delete[] elemSideOrientationsCast;
  delete[] elemBoundariesCast;
  delete[] elemNeighborRanksCast;
  delete[] elemMPIIndicesCast;
  delete[] elemGroupCast;

  // Vertices
  std::map<int, CubeVertex> uniqueVertices;
  transform(vertexMap.begin(),
            vertexMap.end(),
            inserter(uniqueVertices, uniqueVertices.begin()),
            flip_pair<CubeVertex, int>);

  /*  // TODO: not necessary anymore?
    int* vrtxSize = new int[numPartitions[3]];
    std::fill(vrtxSize, vrtxSize + numPartitions[3], uniqueVertices.size());
    checkNcError(nc_put_var_int(ncFile, ncVarVrtxSize, vrtxSize));
    delete[] vrtxSize;
    writes_done++;
    loadBar(writes_done, netcdf_writes);
  */

  m_vertices.resize(uniqueVertices.size());

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

        // checkNcError(nc_put_vara_double(ncFile, ncVarVrtxCoords, start, count, vrtxCoords));
        /*
                for (int i = groupSize - 1; i >= 0; i--) {
                  // MPI_SEND? see NetcdfReader.cpp line 400ff
                }
        */
        writes_done++;
        loadBar(writes_done, netcdf_writes);
      }
    }
  }

  // TODO: communicate the calculation of vrtxCoords, is this if else even necessary?
  if (masterRank >= 0) {
    size_t start[3] = {0, 0, 0};
    size_t count[3] = {1, 0, 3};

    for (int i = groupSize - 1; i >= 0; i--) {
      start[0] = static_cast<size_t>(i + rank);
      count[1] = static_cast<size_t>(sizes[i]);

      // not necessary as we compute this information in the loop above?
//      checkNcError(nc_get_vara_double(
//          ncFile, ncVarVrtxCoords, start, count, reinterpret_cast<double*>(vrtxCoords)));
      if (i != 0) {
#ifdef USE_MPI
        MPI_Send(vrtxCoords, 3 * sizes[i], MPI_DOUBLE, i + rank, 0, seissol::MPI::mpi.comm());
#else  // USE_MPI
      assert(false);
#endif // USE_MPI
      }
    }
  } else {
#ifdef USE_MPI
    MPI_Recv(vrtxCoords,
             3 * sizes[0],
             MPI_DOUBLE,
             master,
             0,
             seissol::MPI::mpi.comm(),
             MPI_STATUS_IGNORE);
#else  // USE_MPI
  assert(false);
#endif // USE_MPI
  }

  

  VrtxCoords* vrtxCoordsCast = reinterpret_cast<VrtxCoords*>(vrtxCoords); 
  // Copy buffers to vertices
  int vrtxSize = sizes[0]; 
  for (int i = 0; i < /*sizes[0]*/uniqueVertices.size(); i++) {
    memcpy(m_vertices[i].coords, &vrtxCoordsCast[i], sizeof(VrtxCoords));
  }

  // delete[] vrtxCoords;
  delete[] vrtxCoordsCast;

  // TODO: Boundaries (MPI neighbours)

  if (masterRank >= 0) {
    for (int i = groupSize - 1; i >= 0; i--) {
      size_t start = static_cast<size_t>(i + rank);

      int size; // TODO: size = bndSize, see nc_put(ncVarBndSize above)
      //checkNcError(nc_get_var1_int(ncFile, ncVarBndSize, &start, &size));
      size = bndSizePtr[start];

      sizes[i] = size;

      if (i != 0) {
#ifdef USE_MPI
        MPI_Send(&size, 1, MPI_INT, i + rank, 0, seissol::MPI::mpi.comm());
#else  // USE_MPI
        assert(false);
#endif // USE_MPI
      }
    }
  } else {
#ifdef USE_MPI
    MPI_Recv(sizes, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
#else  // USE_MPI
    assert(false);
#endif // USE_MPI
  }


  // Get maximum number of neighbors (required to get collective MPI-IO right)
  int maxNeighbors = bndSize;
  // MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, seissol::MPI::mpi.comm());
  int* bndElemLocalIds = new int[bndElemSize];

  //        SCOREP_USER_REGION_DEFINE( r_read_boundaries );
  //        SCOREP_USER_REGION_BEGIN( r_read_boundaries, "read_boundaries",
  //        SCOREP_USER_REGION_TYPE_COMMON );

  size_t bndStart[3] = {0, 0, 0};
  for (int i = 0; i < maxNeighbors; i++) {
    bndStart[1] = static_cast<size_t>(i);

    if (masterRank >= 0) {
      for (int j = groupSize - 1; j >= 0; j--) {
        bndStart[0] = static_cast<size_t>(j + rank);

        // Get neighbor rank from netcdf
        int bndRank;
        //checkNcError(nc_get_var1_int(ncFile, ncVarBndElemRank, bndStart, &bndRank));
        // TODO: was bndStart[0]*numPartitions[3] + bndStart[1]
        bndRank = bndElemRankPtr[bndStart[0]*bndSize + bndStart[1]];

        // Read size of this boundary from netcdf
        int elemSize;
        //checkNcError(nc_get_var1_int(ncFile, ncVarBndElemSize, bndStart, &elemSize)); 
        //TODO: VarBndElemSize has dims [groupSize - 1 + rank][bndSize]
        elemSize = bndElemSizePtr[bndStart[0]*bndSize + bndStart[1]];

        size_t bndCount[3] = {1, 1, bndElemSize};
        //checkNcError(
        //    nc_get_vara_int(ncFile, ncVarBndElemLocalIds, bndStart, bndCount, bndElemLocalIds));
        // TODO: size_t* bndElemLocalIdsPtr[numPartitions[3] * bndSize * bndElemSize]
        memcpy(bndElemLocalIds, &bndElemLocalIdsPtr[(bndStart[0]*bndSize + bndStart[1])*bndElemSize], sizeof(int)*bndCount[2]);

        if (i < sizes[j]) {
          if (j == 0) {
            addMPINeighbor(i, bndRank, elemSize, bndElemLocalIds);
          } else {
#ifdef USE_MPI
            MPI_Send(&bndRank, 1, MPI_INT, j + rank, 0, seissol::MPI::mpi.comm());
            MPI_Send(&elemSize, 1, MPI_INT, j + rank, 0, seissol::MPI::mpi.comm());
            MPI_Send(bndElemLocalIds, elemSize, MPI_INT, j + rank, 0, seissol::MPI::mpi.comm());
#else  // USE_MPI
            assert(false);
#endif // USE_MPI
          }
        }
      }
    } else {
      if (i < sizes[0]) {
#ifdef USE_MPI

        int bndRank;
        MPI_Recv(&bndRank, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
        int elemSize;
        MPI_Recv(&elemSize, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);

        MPI_Recv(bndElemLocalIds,
                 elemSize,
                 MPI_INT,
                 master,
                 0,
                 seissol::MPI::mpi.comm(),
                 MPI_STATUS_IGNORE);

        addMPINeighbor(i, bndRank, elemSize, bndElemLocalIds);
#else  // USE_MPI
        assert(false);
#endif // USE_MPI
      }
    }
  }

  delete[] bndLocalIds;
  delete[] bndElemLocalIds;

  //        SCOREP_USER_REGION_END( r_read_boundaries )

  delete[] sizes;

  // delete all extra pointers
  delete[] bndSizePtr;
  delete[] bndElemSizePtr;
  delete[] bndElemRankPtr;
  delete[] bndElemLocalIdsPtr;
  delete[] elemMPIIndicesPtr;

  // Close netcdf file
  if (masterRank >= 0) {
    checkNcError(nc_close(ncFile));
#ifdef USE_MPI
    MPI_Comm_free(&commMaster);
#endif // USE_MPI
  }

  // Recompute additional information
  findElementsPerVertex();

  logInfo(rank) << "Finished";
}



void seissol::geometry::CubeGenerator::findElementsPerVertex() {
  for (std::vector<Element>::const_iterator i = m_elements.begin(); i != m_elements.end(); i++) {
    for (int j = 0; j < 4; j++) {
      assert(i->vertices[j] < static_cast<int>(m_vertices.size()));
      m_vertices[i->vertices[j]].elements.push_back(i->localId);
    }
  }
}

void seissol::geometry::CubeGenerator::collectiveAccess(int ncFile, int ncVar) {
#ifdef USE_MPI
  checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
#endif // USE_MPI
}

void seissol::geometry::CubeGenerator::addMPINeighbor(int localID,
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
#endif // USE_NETCDF

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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
 */

#include "Geometry/MeshDefinition.h"
#include <Common/Constants.h>
#include <Geometry/MeshReader.h>
#include <Initializer/Parameters/MeshParameters.h>
#include <PUML/TypeInference.h>
#include <PUML/Upward.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <mpi.h>
#include <numeric>
#include <string>

#include "PUMLReader.h"
#include "PartitioningLib.h"

#include "PUML/Downward.h"
#include "PUML/Neighbor.h"
#include "PUML/PUML.h"
#include "PUML/Partition.h"
#include "PUML/PartitionGraph.h"
#include "PUML/PartitionTarget.h"

#include "Monitoring/Instrumentation.h"

#include "Initializer/TimeStepping/LtsWeights/LtsWeights.h"

#include <fstream>
#include <hdf5.h>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace {
/*
 * Possible types of boundary conditions for SeisSol.
 */
enum class BCType { Internal, External, Unknown };

/**
 * Decodes the boundary condition tag into a BCType.
 */
constexpr BCType bcToType(int id) {
  if (id == 0 || id == 3 || id > 64) {
    return BCType::Internal;
  } else if (id == 1 || id == 2 || id == 4 || id == 5 || id == 6 || id == 7) {
    return BCType::External;
  } else {
    return BCType::Unknown;
  }
}

/**
 * Decodes the boundary condition tag into a string representation.
 */
inline std::string bcToString(int id) {
  if (id == 0) {
    return std::string("regular");
  } else if (id == 1) {
    return std::string("free surface");
  } else if (id == 2) {
    return std::string("free surface with gravity");
  } else if (id == 3) {
    return std::string("dynamic rupture");
  } else if (id == 4) {
    return std::string("dirichlet");
  } else if (id == 5) {
    return std::string("absorbing");
  } else if (id == 6) {
    return std::string("periodic");
  } else if (id == 7) {
    return std::string("analytic");
  } else if (id > 64) {
    std::stringstream s;
    s << "fault-tagging (" << id << ")";
    return s.str();
  } else {
    std::stringstream s;
    s << "unknown (" << id << ")";
    return s.str();
  }
}

/**
 * Check if the mesh is locally correct:
 * - if a face is an internal face, the neighbor has to exist.
 * - if a face is an external face, the neighbor must not exist.
 * @param face: face to check
 * @param cellNeighbors: ids of the neighboring cells
 * @param side: side of the tetrahedron to check
 * @param sideBC: boundary condition tag at the side to check
 * @param cellIdAsInFile: Original cell id as it is given in the h5 file
 */
inline bool checkMeshCorrectnessLocally(const PUML::TETPUML::face_t& face,
                                        const int* cellNeighbors,
                                        int side,
                                        int sideBC,
                                        uint64_t cellIdAsInFile) {
  // all of these will only issue warnings here -- the "logError()" is supposed to come later, after
  // all warning have been logged

  // if a face is an internal face, it has to have a neighbor on either this rank or somewhere else:
  if (bcToType(sideBC) == BCType::Internal) {
    if (cellNeighbors[side] < 0 && !face.isShared()) {
      logWarning() << "Element" << cellIdAsInFile << ", side" << side << " has a"
                   << bcToString(sideBC)
                   << "boundary condition, but the neighboring element doesn't exist";
      return false;
    }
  }
  // external boundaries must not have neighboring elements:
  else if (bcToType(sideBC) == BCType::External) {
    if (cellNeighbors[side] >= 0 || face.isShared()) {
      logWarning() << "Element" << cellIdAsInFile << ", side" << side << " has a"
                   << bcToString(sideBC) << "boundary condition, but a neighboring element exists";
      return false;
    }
  }
  // ignore unknown boundary conditions and warn
  else {
    logWarning() << "Element" << cellIdAsInFile << ", side" << side << " has a boundary condition ("
                 << sideBC << ") which is not understood by this version of SeisSol";
    return true;
  }
  return true;
}

// helper arrays

// converts the PUML vertex indexing to the internal SeisSol indexing
const int PumlFaceToSeisSol[4] = {0, 1, 3, 2};

// indexes the vertices on each face i (or FaceVertexToOrientation[i][j] == -1 to indicate that the
// vertex does not lie on it)
const int FaceVertexToOrientation[4][4] = {
    {0, 2, 1, -1}, {0, 1, -1, 2}, {0, -1, 2, 1}, {-1, 0, 1, 2}};

// the first vertex on the face (i.e. FirstFaceVertex[i] == j, where j is the lowest index in
// FaceVertexToOrientation[i] to not be -1)
const int FirstFaceVertex[4] = {0, 0, 0, 1};
} // namespace

seissol::geometry::PUMLReader::PUMLReader(
    const char* meshFile,
    const char* partitioningLib,
    double maximumAllowedTimeStep,
    const char* checkPointFile,
    seissol::initializer::parameters::BoundaryFormat boundaryFormat,
    initializer::time_stepping::LtsWeights* ltsWeights,
    double tpwgt,
    bool readPartitionFromFile)
    : seissol::geometry::MeshReader(MPI::mpi.rank()), boundaryFormat(boundaryFormat) {
  PUML::TETPUML puml;
  puml.setComm(MPI::mpi.comm());

  read(puml, meshFile);

  generatePUML(puml); // We need to call generatePUML in order to create the dual graph of the mesh
  if (ltsWeights != nullptr) {
    ltsWeights->computeWeights(puml, maximumAllowedTimeStep);
  }
  partition(
      puml, ltsWeights, tpwgt, meshFile, partitioningLib, readPartitionFromFile, checkPointFile);

  generatePUML(puml);

  getMesh(puml);
}

void seissol::geometry::PUMLReader::read(PUML::TETPUML& puml, const char* meshFile) {
  SCOREP_USER_REGION("PUMLReader_read", SCOREP_USER_REGION_TYPE_FUNCTION);

  const std::string file(meshFile);

  puml.open((file + ":/connect").c_str(), (file + ":/geometry").c_str());
  puml.addData<int>((file + ":/group").c_str(), PUML::CELL, {});

  if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I32) {
    puml.addData<uint32_t>((file + ":/boundary").c_str(), PUML::CELL, {});
  } else if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I64) {
    puml.addData<uint64_t>((file + ":/boundary").c_str(), PUML::CELL, {});
  } else if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I32x4) {
    puml.addData<int>((file + ":/boundary").c_str(), PUML::CELL, {4});
  }

  // TODO(David): change to uint64_t/size_t once we have an MPI module for that ready
  const size_t localCells = puml.numOriginalCells();
  size_t localStart = 0;

#ifdef USE_MPI
  MPI_Exscan(&localCells, &localStart, 1, PUML::MPITypeInfer<size_t>::type(), MPI_SUM, puml.comm());
#endif

  std::vector<size_t> cellIdsAsInFile(localCells);
  std::iota(cellIdsAsInFile.begin(), cellIdsAsInFile.end(), localStart);
  puml.addDataArray(cellIdsAsInFile.data(), PUML::CELL, {});
}

int seissol::geometry::PUMLReader::readPartition(PUML::TETPUML& puml,
                                                 int* partition,
                                                 const char* checkPointFile) {
  /*
  write the partionning array to an hdf5 file using parallel access
  see https://support.hdfgroup.org/ftp/HDF5/examples/parallel/coll_test.c for more info about the
  hdf5 functions
  */
  SCOREP_USER_REGION("PUMLReader_readPartition", SCOREP_USER_REGION_TYPE_FUNCTION);
  const int rank = seissol::MPI::mpi.rank();
  const int nrank = seissol::MPI::mpi.size();
  int nPartitionCells = puml.numOriginalCells();

  /*
   Gather number of cells in each nodes. This is necessary to be able to write the data in the
   correct location
  */
  int* numCells = new int[nrank];
  int* offsets = new int[nrank];
  MPI_Allgather(&nPartitionCells, 1, MPI_INT, numCells, 1, MPI_INT, MPI::mpi.comm());

  offsets[0] = 0;
  for (int rk = 1; rk < nrank; ++rk) {
    offsets[rk] = offsets[rk - 1] + numCells[rk - 1];
  }
  const hsize_t dimMem[] = {static_cast<hsize_t>(nPartitionCells)};

  /*
   Open file and dataset
  */
  MPI_Info info = MPI_INFO_NULL;
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plistId, seissol::MPI::mpi.comm(), info);

  std::ostringstream os;
  os << checkPointFile << "_partitions_o" << ConvergenceOrder << "_n" << nrank << ".h5";
  const std::string fname = os.str();

  const std::ifstream ifile(fname.c_str());
  if (!ifile) {
    logInfo(rank) << fname.c_str() << "does not exist";
    return -1;
  }

  const hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plistId);
  H5Pclose(plistId);

  const hid_t dataset = H5Dopen2(file, "/partition", H5P_DEFAULT);
  /*
   Create memspace (portion of filespace) and read collectively the data
  */
  const hid_t memspace = H5Screate_simple(1, dimMem, nullptr);
  const hid_t filespace = H5Dget_space(dataset);

  hsize_t start[] = {static_cast<hsize_t>(offsets[rank])};
  hsize_t count[] = {static_cast<hsize_t>(nPartitionCells)};
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  plistId = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

  const int status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace, plistId, partition);

  if (status < 0)
    logError() << "An error occured when reading the partitionning with HDF5";
  H5Dclose(dataset);
  H5Fclose(file);

  logInfo(rank) << "partitionning was read successfully from " << fname.c_str();
  return 0;
}

void seissol::geometry::PUMLReader::writePartition(PUML::TETPUML& puml,
                                                   int* partition,
                                                   const char* checkPointFile) {
  /*
  write the partionning array to an hdf5 file using parallel access
  see https://support.hdfgroup.org/ftp/HDF5/examples/parallel/coll_test.c for more info about the
  hdf5 functions
  */
  SCOREP_USER_REGION("PUMLReader_writePartition", SCOREP_USER_REGION_TYPE_FUNCTION);
  const int rank = seissol::MPI::mpi.rank();
  const int nrank = seissol::MPI::mpi.size();
  int nPartitionCells = puml.numOriginalCells();

  /*
   Gather number of cells in each nodes. This is necessary to be able to write the data in the
   correct location
  */
  int* numCells = new int[nrank];
  int* offsets = new int[nrank];
  MPI_Allgather(&nPartitionCells, 1, MPI_INT, numCells, 1, MPI_INT, MPI::mpi.comm());

  offsets[0] = 0;
  for (int rk = 1; rk < nrank; ++rk) {
    offsets[rk] = offsets[rk - 1] + numCells[rk - 1];
  }
  const int nCells = offsets[nrank - 1] + numCells[nrank - 1];

  const hsize_t dim[] = {static_cast<hsize_t>(nCells)};
  const hsize_t dimMem[] = {static_cast<hsize_t>(nPartitionCells)};

  /*
   Create file and file space
  */
  MPI_Info info = MPI_INFO_NULL;
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plistId, seissol::MPI::mpi.comm(), info);

  std::ostringstream os;
  os << checkPointFile << "_partitions_o" << ConvergenceOrder << "_n" << nrank << ".h5";
  const std::string fname = os.str();

  const hid_t file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
  H5Pclose(plistId);

  hid_t filespace = H5Screate_simple(1, dim, nullptr);
  const hid_t dataset = H5Dcreate(
      file, "/partition", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(filespace);

  /*
   Create memspace (portion of filespace) and write collectively the data
  */
  const hid_t memspace = H5Screate_simple(1, dimMem, nullptr);
  filespace = H5Dget_space(dataset);

  hsize_t start[] = {static_cast<hsize_t>(offsets[rank])};
  hsize_t count[] = {static_cast<hsize_t>(nPartitionCells)};
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  plistId = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plistId, H5FD_MPIO_COLLECTIVE);

  const int status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plistId, partition);

  if (status < 0)
    logError() << "An error occured when writing the partitionning with HDF5";
  H5Dclose(dataset);
  H5Fclose(file);
}

void seissol::geometry::PUMLReader::partition(PUML::TETPUML& puml,
                                              initializer::time_stepping::LtsWeights* ltsWeights,
                                              double tpwgt,
                                              const char* meshFile,
                                              const char* partitioningLib,
                                              bool readPartitionFromFile,
                                              const char* checkPointFile) {
  SCOREP_USER_REGION("PUMLReader_partition", SCOREP_USER_REGION_TYPE_FUNCTION);

  auto doPartition =
      [&] {
        auto partType = toPartitionerType(std::string_view(partitioningLib));
        logInfo(MPI::mpi.rank()) << "Using the" << toStringView(partType)
                                 << "partition library and strategy.";
        if (partType == PUML::PartitionerType::None) {
          logWarning(MPI::mpi.rank())
              << partitioningLib
              << "not found. Expect poor performance as the mesh is not properly partitioned.";
        }
        auto partitioner = PUML::TETPartition::getPartitioner(partType);
        if (partitioner == nullptr) {
          logError() << "Unrecognized partition library: " << partitioningLib;
        }
        auto graph = PUML::TETPartitionGraph(puml);
        graph.setVertexWeights(ltsWeights->vertexWeights(), ltsWeights->nWeightsPerVertex());

#ifdef USE_MPI
        auto nodeWeights = std::vector<double>(MPI::mpi.size());
        MPI_Allgather(
            &tpwgt, 1, MPI_DOUBLE, nodeWeights.data(), 1, MPI_DOUBLE, seissol::MPI::mpi.comm());
        double sum = 0.0;
        for (const auto& w : nodeWeights) {
          sum += w;
        }
        for (auto& w : nodeWeights) {
          w /= sum;
        }
#else
        auto nodeWeights = std::vector<double>{1.0};
#endif

        auto target = PUML::PartitionTarget{};
        target.setVertexWeights(nodeWeights);
        target.setImbalance(ltsWeights->imbalances()[0] - 1.0);

        return partitioner->partition(graph, target);
      };

  auto newPartition = std::vector<int>();
  if (readPartitionFromFile) {
    newPartition.resize(puml.numOriginalCells());
    const int status = readPartition(puml, newPartition.data(), checkPointFile);
    if (status < 0) {
      newPartition = doPartition();
      writePartition(puml, newPartition.data(), checkPointFile);
    }
  } else {
    newPartition = doPartition();
  }

  puml.partition(newPartition.data());
}

void seissol::geometry::PUMLReader::generatePUML(PUML::TETPUML& puml) {
  SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

  puml.generateMesh();
}

void seissol::geometry::PUMLReader::getMesh(const PUML::TETPUML& puml) {
  SCOREP_USER_REGION("PUMLReader_getmesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const int rank = MPI::mpi.rank();

  const std::vector<PUML::TETPUML::cell_t>& cells = puml.cells();
  const std::vector<PUML::TETPUML::face_t>& faces = puml.faces();
  const std::vector<PUML::TETPUML::vertex_t>& vertices = puml.vertices();

  const int* group = reinterpret_cast<const int*>(puml.cellData(0));
  const void* boundaryCond = puml.cellData(1);
  const auto* cellIdsAsInFile = reinterpret_cast<const size_t*>(puml.cellData(2));

  std::unordered_map<int, std::vector<unsigned int>> neighborInfo; // List of shared local face ids

  bool isMeshCorrect = true;

  // Compute everything local
  m_elements.resize(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    m_elements[i].globalId = cellIdsAsInFile[i];
    m_elements[i].localId = i;

    // Vertices
    PUML::Downward::vertices(
        puml, cells[i], reinterpret_cast<unsigned int*>(m_elements[i].vertices));

    // Neighbor information
    unsigned int faceids[4];
    PUML::Downward::faces(puml, cells[i], faceids);
    int neighbors[4];
    PUML::Neighbor::face(puml, i, neighbors);
    for (int j = 0; j < 4; j++) {
      int bcCurrentFace = decodeBoundary(boundaryCond, i, j, boundaryFormat);
      const bool isLocallyCorrect = checkMeshCorrectnessLocally(
          faces[faceids[j]], neighbors, j, bcCurrentFace, cellIdsAsInFile[i]);
      isMeshCorrect &= isLocallyCorrect;
      if (neighbors[j] < 0) {
        m_elements[i].neighbors[PumlFaceToSeisSol[j]] = cells.size();

        if (!faces[faceids[j]].isShared()) {
          // Boundary sides
          m_elements[i].neighborRanks[PumlFaceToSeisSol[j]] = rank;
        } else {
          // MPI Boundary
          neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

          m_elements[i].neighborRanks[PumlFaceToSeisSol[j]] = faces[faceids[j]].shared()[0];
        }
      } else {
        assert(neighbors[j] >= 0 && static_cast<unsigned>(neighbors[j]) < cells.size());

        m_elements[i].neighbors[PumlFaceToSeisSol[j]] = neighbors[j];

        int nfaces[4];
        PUML::Neighbor::face(puml, neighbors[j], nfaces);
        int* back = std::find(nfaces, nfaces + 4, i);
        assert(back < nfaces + 4);

        m_elements[i].neighborSides[PumlFaceToSeisSol[j]] = PumlFaceToSeisSol[back - nfaces];

        const unsigned int firstVertex =
            m_elements[i].vertices[FirstFaceVertex[PumlFaceToSeisSol[j]]];

        unsigned int nvertices[4];
        PUML::Downward::vertices(puml, cells[neighbors[j]], nvertices);
        unsigned int* neighborFirstVertex = std::find(nvertices, nvertices + 4, firstVertex);

        m_elements[i].sideOrientations[PumlFaceToSeisSol[j]] =
            FaceVertexToOrientation[m_elements[i].neighborSides[PumlFaceToSeisSol[j]]]
                                   [neighborFirstVertex - nvertices];
        assert(m_elements[i].sideOrientations[PumlFaceToSeisSol[j]] >= 0);

        m_elements[i].neighborRanks[PumlFaceToSeisSol[j]] = rank;
      }

      const int faultTag = bcCurrentFace;
      if (bcCurrentFace > 64) {
        bcCurrentFace = 3;
      }
      m_elements[i].boundaries[PumlFaceToSeisSol[j]] = bcCurrentFace;
      m_elements[i].faultTags[PumlFaceToSeisSol[j]] = faultTag;
      m_elements[i].mpiIndices[PumlFaceToSeisSol[j]] = 0;
    }

    m_elements[i].group = group[i];
  }
  if (!isMeshCorrect) {
    logError() << "Found at least one broken face in the mesh, see errors above for a more "
                  "detailled analysis.";
  }

  // Exchange ghost layer information and generate neighbor list
  char** copySide = new char*[neighborInfo.size()];
  char** ghostSide = new char*[neighborInfo.size()];
  auto** copyFirstVertex = new unsigned long*[neighborInfo.size()];
  auto** ghostFirstVertex = new unsigned long*[neighborInfo.size()];

  auto* requests = new MPI_Request[neighborInfo.size() * 4];

  std::unordered_set<unsigned int> t;
#ifndef NDEBUG
  unsigned int sum = 0;
#endif
  unsigned int k = 0;
  for (auto it = neighborInfo.begin(); it != neighborInfo.end(); ++it, ++k) {
    // Need to sort the neighborInfo vectors once
    std::sort(it->second.begin(), it->second.end(), [&](unsigned int a, unsigned int b) {
      return puml.faces()[a].gid() < puml.faces()[b].gid();
    });

    t.insert(it->second.begin(), it->second.end());
#ifndef NDEBUG
    sum += it->second.size();
#endif

    // Create MPI neighbor list
    addMPINeighor(puml, it->first, it->second);

    copySide[k] = new char[it->second.size()];
    ghostSide[k] = new char[it->second.size()];
    copyFirstVertex[k] = new unsigned long[it->second.size()];
    ghostFirstVertex[k] = new unsigned long[it->second.size()];

    MPI_Irecv(
        ghostSide[k], it->second.size(), MPI_CHAR, it->first, 0, MPI::mpi.comm(), &requests[k]);
    MPI_Irecv(ghostFirstVertex[k],
              it->second.size(),
              MPI_UNSIGNED_LONG,
              it->first,
              0,
              MPI::mpi.comm(),
              &requests[neighborInfo.size() + k]);

    // Neighbor side
    for (unsigned int i = 0; i < it->second.size(); i++) {
      // The side of boundary
      int cellIds[2];
      PUML::Upward::cells(puml, faces[it->second[i]], cellIds);
      const int side = PUML::Downward::faceSide(puml, cells[cellIds[0]], it->second[i]);
      assert(side >= 0 && side < 4);
      copySide[k][i] = side;

      // First vertex of the face on the boundary
      const unsigned int firstVertex =
          m_elements[cellIds[0]].vertices[FirstFaceVertex[PumlFaceToSeisSol[side]]];
      copyFirstVertex[k][i] = vertices[firstVertex].gid();

      // Set the MPI index
      assert(m_elements[cellIds[0]].mpiIndices[PumlFaceToSeisSol[side]] == 0);
      m_elements[cellIds[0]].mpiIndices[PumlFaceToSeisSol[side]] = i;
    }

    MPI_Isend(copySide[k],
              it->second.size(),
              MPI_CHAR,
              it->first,
              0,
              MPI::mpi.comm(),
              &requests[neighborInfo.size() * 2 + k]);
    MPI_Isend(copyFirstVertex[k],
              it->second.size(),
              MPI_UNSIGNED_LONG,
              it->first,
              0,
              MPI::mpi.comm(),
              &requests[neighborInfo.size() * 3 + k]);
  }
#ifndef NDEBUG
  assert(t.size() == sum);
#endif

  MPI_Waitall(neighborInfo.size() * 4, requests, MPI_STATUSES_IGNORE);

  k = 0;
  for (auto it = neighborInfo.begin(); it != neighborInfo.end(); ++it, k++) {
    for (unsigned int i = 0; i < it->second.size(); i++) {
      // Set neighbor side
      int cellIds[2];
      PUML::Upward::cells(puml, faces[it->second[i]], cellIds);
      assert(cellIds[1] < 0);

      const int side = copySide[k][i];
      const int gSide = ghostSide[k][i];
      m_elements[cellIds[0]].neighborSides[PumlFaceToSeisSol[side]] = PumlFaceToSeisSol[gSide];

      // Set side sideOrientation
      unsigned long nvertices[4];
      PUML::Downward::gvertices(puml, cells[cellIds[0]], nvertices);

      unsigned long* localFirstVertex = std::find(nvertices, nvertices + 4, ghostFirstVertex[k][i]);
      assert(localFirstVertex != nvertices + 4);

      m_elements[cellIds[0]].sideOrientations[PumlFaceToSeisSol[side]] =
          FaceVertexToOrientation[PumlFaceToSeisSol[side]][localFirstVertex - nvertices];
      assert(m_elements[cellIds[0]].sideOrientations[PumlFaceToSeisSol[side]] >= 0);
    }

    delete[] copySide[k];
    delete[] ghostSide[k];
    delete[] copyFirstVertex[k];
    delete[] ghostFirstVertex[k];
  }

  delete[] copySide;
  delete[] ghostSide;
  delete[] copyFirstVertex;
  delete[] ghostFirstVertex;
  delete[] requests;

  // Set vertices
  m_vertices.resize(vertices.size());
  for (std::size_t i = 0; i < vertices.size(); i++) {
    memcpy(m_vertices[i].coords, vertices[i].coordinate(), 3 * sizeof(double));

    const std::vector<int> elementsInt;

    PUML::Upward::cells(puml, vertices[i], m_vertices[i].elements);
  }
}

void seissol::geometry::PUMLReader::addMPINeighor(const PUML::TETPUML& puml,
                                                  int rank,
                                                  const std::vector<unsigned int>& faces) {
  const std::size_t id = m_MPINeighbors.size();
  MPINeighbor& neighbor = m_MPINeighbors[rank];

  neighbor.localID = id;
  neighbor.elements.resize(faces.size());

  for (std::size_t i = 0; i < faces.size(); i++) {
    int cellIds[2];
    PUML::Upward::cells(puml, puml.faces()[faces[i]], cellIds);

    neighbor.elements[i].localElement = cellIds[0];
  }
}

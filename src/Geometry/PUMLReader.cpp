// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "Geometry/MeshDefinition.h"
#include <Common/Iterator.h>
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

namespace seissol::geometry {

PUMLReader::PUMLReader(const char* meshFile,
                       const char* partitioningLib,
                       seissol::initializer::parameters::BoundaryFormat boundaryFormat,
                       initializer::time_stepping::LtsWeights* ltsWeights,
                       double tpwgt)
    : MeshReader(MPI::mpi.rank()), boundaryFormat(boundaryFormat) {
  PUML::TETPUML puml;
  puml.setComm(MPI::mpi.comm());

  read(puml, meshFile);

  generatePUML(puml); // We need to call generatePUML in order to create the dual graph of the mesh
  if (ltsWeights != nullptr) {
    ltsWeights->computeWeights(puml);
  }
  partition(puml, ltsWeights, tpwgt, meshFile, partitioningLib);

  generatePUML(puml);

  getMesh(puml);
}

void PUMLReader::read(PUML::TETPUML& puml, const char* meshFile) {
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

  const size_t localCells = puml.numOriginalCells();
  size_t localStart = 0;

#ifdef USE_MPI
  MPI_Exscan(&localCells, &localStart, 1, PUML::MPITypeInfer<size_t>::type(), MPI_SUM, puml.comm());
#endif

  std::vector<size_t> cellIdsAsInFile(localCells);
  std::iota(cellIdsAsInFile.begin(), cellIdsAsInFile.end(), localStart);
  puml.addDataArray(cellIdsAsInFile.data(), PUML::CELL, {});
}

void PUMLReader::partition(PUML::TETPUML& puml,
                           initializer::time_stepping::LtsWeights* ltsWeights,
                           double tpwgt,
                           const char* meshFile,
                           const char* partitioningLib) {
  SCOREP_USER_REGION("PUMLReader_partition", SCOREP_USER_REGION_TYPE_FUNCTION);

  auto partType = toPartitionerType(std::string_view(partitioningLib));
  logInfo() << "Using the" << toStringView(partType) << "partition library and strategy.";
  if (partType == PUML::PartitionerType::None) {
    logWarning() << partitioningLib
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
  MPI_Allgather(&tpwgt, 1, MPI_DOUBLE, nodeWeights.data(), 1, MPI_DOUBLE, seissol::MPI::mpi.comm());
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

  auto newPartition = partitioner->partition(graph, target);

  puml.addDataArray(ltsWeights->clusterIds().data(), PUML::CELL, {});
  puml.addDataArray(ltsWeights->timesteps().data(), PUML::CELL, {});

  puml.partition(newPartition.data());
}

void PUMLReader::generatePUML(PUML::TETPUML& puml) {
  SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

  puml.generateMesh();
}

void PUMLReader::getMesh(const PUML::TETPUML& puml) {
  SCOREP_USER_REGION("PUMLReader_getmesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const int rank = MPI::mpi.rank();

  const std::vector<PUML::TETPUML::cell_t>& cells = puml.cells();
  const std::vector<PUML::TETPUML::face_t>& faces = puml.faces();
  const std::vector<PUML::TETPUML::vertex_t>& vertices = puml.vertices();

  const int* group = reinterpret_cast<const int*>(puml.cellData(0));
  const void* boundaryCond = puml.cellData(1);
  const auto* cellIdsAsInFile = reinterpret_cast<const size_t*>(puml.cellData(2));
  const auto* clusterIds = reinterpret_cast<const int*>(puml.cellData(3));
  const auto* timestep = reinterpret_cast<const double*>(puml.cellData(4));

  std::unordered_map<int, std::vector<unsigned int>> neighborInfo; // List of shared local face ids

  bool isMeshCorrect = true;

  // Compute everything local
  m_elements.resize(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    m_elements[i].globalId = cellIdsAsInFile[i];
    m_elements[i].localId = i;
    m_elements[i].clusterId = clusterIds[i];
    m_elements[i].timestep = timestep[i];

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
  std::vector<std::vector<char>> copySide(neighborInfo.size());
  std::vector<std::vector<char>> ghostSide(neighborInfo.size());
  std::vector<std::vector<unsigned long>> copyFirstVertex(neighborInfo.size());
  std::vector<std::vector<unsigned long>> ghostFirstVertex(neighborInfo.size());

  std::vector<MPI_Request> requests(neighborInfo.size() * 4);

  std::unordered_set<unsigned int> t;
#ifndef NDEBUG
  unsigned int sum = 0;
#endif
  for (auto [k, info] : seissol::common::enumerate(neighborInfo)) {
    // Need to sort the neighborInfo vectors once
    std::sort(info.second.begin(), info.second.end(), [&](unsigned int a, unsigned int b) {
      return puml.faces()[a].gid() < puml.faces()[b].gid();
    });

    t.insert(info.second.begin(), info.second.end());
#ifndef NDEBUG
    sum += info.second.size();
#endif

    // Create MPI neighbor list
    addMPINeighor(puml, info.first, info.second);

    copySide[k].resize(info.second.size());
    ghostSide[k].resize(info.second.size());
    copyFirstVertex[k].resize(info.second.size());
    ghostFirstVertex[k].resize(info.second.size());

    MPI_Irecv(ghostSide[k].data(),
              info.second.size(),
              MPI_CHAR,
              info.first,
              0,
              MPI::mpi.comm(),
              &requests[k]);
    MPI_Irecv(ghostFirstVertex[k].data(),
              info.second.size(),
              MPI_UNSIGNED_LONG,
              info.first,
              0,
              MPI::mpi.comm(),
              &requests[neighborInfo.size() + k]);

    // Neighbor side
    for (unsigned int i = 0; i < info.second.size(); i++) {
      // The side of boundary
      int cellIds[2];
      PUML::Upward::cells(puml, faces[info.second[i]], cellIds);
      const int side = PUML::Downward::faceSide(puml, cells[cellIds[0]], info.second[i]);
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

    MPI_Isend(copySide[k].data(),
              info.second.size(),
              MPI_CHAR,
              info.first,
              0,
              MPI::mpi.comm(),
              &requests[neighborInfo.size() * 2 + k]);
    MPI_Isend(copyFirstVertex[k].data(),
              info.second.size(),
              MPI_UNSIGNED_LONG,
              info.first,
              0,
              MPI::mpi.comm(),
              &requests[neighborInfo.size() * 3 + k]);
  }
#ifndef NDEBUG
  assert(t.size() == sum);
#endif

  MPI_Waitall(neighborInfo.size() * 4, requests.data(), MPI_STATUSES_IGNORE);

  for (auto [k, info] : seissol::common::enumerate(neighborInfo)) {
    for (unsigned int i = 0; i < info.second.size(); i++) {
      // Set neighbor side
      int cellIds[2];
      PUML::Upward::cells(puml, faces[info.second[i]], cellIds);
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
  }

  // Set vertices
  m_vertices.resize(vertices.size());
  for (std::size_t i = 0; i < vertices.size(); i++) {
    memcpy(m_vertices[i].coords, vertices[i].coordinate(), 3 * sizeof(double));

    PUML::Upward::cells(puml, vertices[i], m_vertices[i].elements);
  }
}

void PUMLReader::addMPINeighor(const PUML::TETPUML& puml,
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

bool PUMLReader::inlineTimestepCompute() const { return true; }

bool PUMLReader::inlineClusterCompute() const { return true; }

} // namespace seissol::geometry

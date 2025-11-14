// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#include "PUMLReader.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Geometry/MeshDefinition.h"
#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Initializer/TimeStepping/LtsWeights/LtsWeights.h"
#include "Monitoring/Instrumentation.h"
#include "PartitioningLib.h"

#include <PUML/Downward.h>
#include <PUML/Neighbor.h>
#include <PUML/PUML.h>
#include <PUML/Partition.h>
#include <PUML/PartitionGraph.h>
#include <PUML/PartitionTarget.h>
#include <PUML/Topology.h>
#include <PUML/TypeInference.h>
#include <PUML/Upward.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <hdf5.h>
#include <mpi.h>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace {

// PUML sanity checks
using PumlTopologyInternal = PUML::internal::Topology<seissol::geometry::PumlTopology>;
static_assert(PumlTopologyInternal::dimension() == seissol::Cell::Dim);
static_assert(PumlTopologyInternal::cellfaces() == seissol::Cell::NumFaces);
static_assert(PumlTopologyInternal::cellvertices() == seissol::Cell::NumVertices);

void logassertI(bool condition, const std::string& file, int line) {
  if (!condition) {
    logError() << "Assertion failure in" << file << "at" << line;
  }
}

#define logassert(x) logassertI(x, __FILE__, __LINE__)

/*
 * Possible types of boundary conditions for SeisSol.
 */
enum class BCType { Internal, External, Unknown };

/**
 * Decodes the boundary condition tag into a BCType.
 */
constexpr BCType bcToType(int id) {
  if (id == 0 || id == 3 || id == 6 || id > 64) {
    return BCType::Internal;
  } else if (id == 1 || id == 2 || id == 4 || id == 5 || id == 7) {
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
template <PUML::TopoType Topo>
inline bool
    checkMeshCorrectnessLocally(const typename PUML::PUML<Topo>::face_t& face,
                                const std::array<int, seissol::Cell::NumFaces>& cellNeighbors,
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
const std::array<std::int32_t, 4> PumlFaceToSeisSol = {0, 1, 3, 2};

// indexes the vertices on each face i (or FaceVertexToOrientation[i][j] == -1 to indicate that the
// vertex does not lie on it)
const std::array<std::array<std::int32_t, 4>, 4> FaceVertexToOrientation = {
    std::array<std::int32_t, 4>{0, 2, 1, -1},
    std::array<std::int32_t, 4>{0, 1, -1, 2},
    std::array<std::int32_t, 4>{0, -1, 2, 1},
    std::array<std::int32_t, 4>{-1, 0, 1, 2}};

// the first vertex on the face (i.e. FirstFaceVertex[i] == j, where j is the lowest index in
// FaceVertexToOrientation[i] to not be -1)
const std::array<std::int32_t, 4> FirstFaceVertex = {0, 0, 0, 1};
} // namespace

namespace seissol::geometry {

PUMLReader::PUMLReader(const std::string& meshFile,
                       const std::string& partitioningLib,
                       seissol::initializer::parameters::BoundaryFormat boundaryFormat,
                       seissol::initializer::parameters::TopologyFormat topologyFormat,
                       initializer::time_stepping::LtsWeights* ltsWeights,
                       double tpwgt)
    : MeshReader(seissol::Mpi::mpi.rank()), boundaryFormat(boundaryFormat) {
  // we need up to two meshes, potentially:
  // one mesh for the geometry
  // one mesh for the topology
  // they will only differ if we have periodic boundary conditions

  PumlMesh meshTopologyExtra;
  PumlMesh meshGeometry;
  meshTopologyExtra.setComm(seissol::Mpi::mpi.comm());
  meshGeometry.setComm(seissol::Mpi::mpi.comm());

  read(meshGeometry, meshFile, false);

  // Note: we need to call generatePUML in order to create the dual graph of the mesh
  // Note 2: we also need it for vertex identification
  meshGeometry.generateMesh();

  if (topologyFormat != initializer::parameters::TopologyFormat::Geometric) {
    // we have a topology mesh; separate from the physical mesh

    read(meshTopologyExtra,
         meshFile,
         topologyFormat == initializer::parameters::TopologyFormat::IdentifyFace);

    int id = -1;
    if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyVertex) {
      id = meshTopologyExtra.addData<unsigned long>(
          (std::string(meshFile) + ":/identify").c_str(), PUML::VERTEX, {});
    }

    // generate the topology mesh for the dual graph
    meshTopologyExtra.generateMesh();

    if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyVertex) {
      // re-identify vertices; then re-distribute
      meshTopologyExtra.identify(id);
      meshTopologyExtra.generateMesh();
    }
  }

  auto& meshTopology = topologyFormat == initializer::parameters::TopologyFormat::Geometric
                           ? meshGeometry
                           : meshTopologyExtra;

  if (ltsWeights != nullptr) {
    ltsWeights->computeWeights(meshTopology, meshGeometry);
  }
  partition(meshTopology, meshGeometry, ltsWeights, tpwgt, partitioningLib);

  generatePUML(meshTopology, meshGeometry);

  getMesh(meshTopology, meshGeometry);
}

void PUMLReader::read(PumlMesh& meshTopology, const std::string& file, bool topology) {
  SCOREP_USER_REGION("PUMLReader_read", SCOREP_USER_REGION_TYPE_FUNCTION);

  if (topology) {
    meshTopology.open((file + ":/topology").c_str(), (file + ":/geometry").c_str());
  } else {
    meshTopology.open((file + ":/connect").c_str(), (file + ":/geometry").c_str());
  }
  meshTopology.addData<int>((file + ":/group").c_str(), PUML::CELL, {});

  if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I32) {
    meshTopology.addData<uint32_t>((file + ":/boundary").c_str(), PUML::CELL, {});
  } else if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I64) {
    meshTopology.addData<uint64_t>((file + ":/boundary").c_str(), PUML::CELL, {});
  } else if (boundaryFormat == seissol::initializer::parameters::BoundaryFormat::I32x4) {
    meshTopology.addData<int>((file + ":/boundary").c_str(), PUML::CELL, {4});
  }

  const size_t localCells = meshTopology.numOriginalCells();
  size_t localStart = 0;

  MPI_Exscan(&localCells,
             &localStart,
             1,
             PUML::MPITypeInfer<size_t>::type(),
             MPI_SUM,
             meshTopology.comm());

  std::vector<size_t> cellIdsAsInFile(localCells);
  std::iota(cellIdsAsInFile.begin(), cellIdsAsInFile.end(), localStart);
  meshTopology.addDataArray(cellIdsAsInFile.data(), PUML::CELL, {});
}

void PUMLReader::partition(PumlMesh& meshTopology,
                           PumlMesh& meshGeometry,
                           initializer::time_stepping::LtsWeights* ltsWeights,
                           double tpwgt,
                           const std::string& partitioningLib) {
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
  auto graph = PUML::TETPartitionGraph(meshTopology);
  graph.setVertexWeights(ltsWeights->vertexWeights(), ltsWeights->nWeightsPerVertex());

  auto nodeWeights = std::vector<double>(Mpi::mpi.size());
  MPI_Allgather(&tpwgt, 1, MPI_DOUBLE, nodeWeights.data(), 1, MPI_DOUBLE, seissol::Mpi::mpi.comm());
  double sum = 0.0;
  for (const auto& w : nodeWeights) {
    sum += w;
  }
  for (auto& w : nodeWeights) {
    w /= sum;
  }

  auto target = PUML::PartitionTarget{};
  target.setVertexWeights(nodeWeights);
  target.setImbalance(ltsWeights->imbalances()[0] - 1.0);

  auto newPartition = partitioner->partition(graph, target);

  meshGeometry.addDataArray(ltsWeights->clusterIds().data(), PUML::CELL, {});
  meshGeometry.addDataArray(ltsWeights->timesteps().data(), PUML::CELL, {});

  meshGeometry.partition(newPartition.data());
  if (&meshTopology != &meshGeometry) {
    meshTopology.partition(newPartition.data());
  }
}

void PUMLReader::generatePUML(PumlMesh& meshTopology, PumlMesh& meshGeometry) {
  SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

  if (&meshTopology != &meshGeometry) {
    meshTopology.generateMesh();
  }
  meshGeometry.generateMesh();
}

void PUMLReader::getMesh(const PumlMesh& meshTopology, const PumlMesh& meshGeometry) {
  SCOREP_USER_REGION("PUMLReader_getmesh", SCOREP_USER_REGION_TYPE_FUNCTION);

  const int rank = Mpi::mpi.rank();

  const std::vector<PumlMesh::cell_t>& cells = meshTopology.cells();
  const std::vector<PumlMesh::face_t>& faces = meshTopology.faces();
  const std::vector<PumlMesh::vertex_t>& vertices = meshTopology.vertices();

  const std::vector<PumlMesh::cell_t>& cellsGeometry = meshGeometry.cells();
  const std::vector<PumlMesh::vertex_t>& verticesGeometry = meshGeometry.vertices();

  const int* group = reinterpret_cast<const int*>(meshGeometry.cellData(0));
  const void* boundaryCond = meshGeometry.cellData(1);
  const auto* cellIdsAsInFile = reinterpret_cast<const size_t*>(meshGeometry.cellData(2));
  const auto* clusterIds = reinterpret_cast<const int*>(meshGeometry.cellData(3));
  const auto* timestep = reinterpret_cast<const double*>(meshGeometry.cellData(4));

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
        meshGeometry, cellsGeometry[i], reinterpret_cast<unsigned int*>(m_elements[i].vertices));

    std::array<unsigned int, Cell::NumVertices> topoVertices{};
    PUML::Downward::vertices(meshTopology, cells[i], topoVertices.data());

    // Neighbor information
    std::array<unsigned int, Cell::NumFaces> faceids{};
    PUML::Downward::faces(meshTopology, cells[i], faceids.data());
    std::array<int, Cell::NumFaces> neighbors{};
    PUML::Neighbor::face(meshTopology, i, neighbors.data());

    for (std::size_t j = 0; j < Cell::NumFaces; j++) {
      int bcCurrentFace = decodeBoundary(boundaryCond, i, j, boundaryFormat);
      const bool isLocallyCorrect = checkMeshCorrectnessLocally<PumlTopology>(
          faces[faceids[j]], neighbors, j, bcCurrentFace, cellIdsAsInFile[i]);
      isMeshCorrect &= isLocallyCorrect;
      if (neighbors[j] < 0) {
        m_elements[i].neighbors[PumlFaceToSeisSol[j]] = cellsGeometry.size();

        if (!faces[faceids[j]].isShared()) {
          // Boundary sides
          m_elements[i].neighborRanks[PumlFaceToSeisSol[j]] = rank;
        } else {
          // MPI Boundary
          neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

          m_elements[i].neighborRanks[PumlFaceToSeisSol[j]] = faces[faceids[j]].shared()[0];
        }
      } else {
        logassert(neighbors[j] >= 0 && static_cast<std::size_t>(neighbors[j]) < cells.size());

        m_elements[i].neighbors[PumlFaceToSeisSol[j]] = neighbors[j];

        std::array<int, Cell::NumFaces> nfaces{};
        PUML::Neighbor::face(meshTopology, neighbors[j], nfaces.data());
        const auto* back = std::find(nfaces.begin(), nfaces.end(), i);
        logassert(back != nfaces.end());

        m_elements[i].neighborSides[PumlFaceToSeisSol[j]] =
            PumlFaceToSeisSol[back - nfaces.begin()];

        const auto firstVertex = topoVertices[FirstFaceVertex[PumlFaceToSeisSol[j]]];

        std::array<unsigned int, Cell::NumVertices> nvertices{};
        PUML::Downward::vertices(meshTopology, cells[neighbors[j]], nvertices.data());
        const auto* neighborFirstVertex =
            std::find(nvertices.begin(), nvertices.end(), firstVertex);
        logassert(neighborFirstVertex != nvertices.end());

        m_elements[i].sideOrientations[PumlFaceToSeisSol[j]] =
            FaceVertexToOrientation[m_elements[i].neighborSides[PumlFaceToSeisSol[j]]]
                                   [neighborFirstVertex - nvertices.begin()];
        logassert(m_elements[i].sideOrientations[PumlFaceToSeisSol[j]] >= 0);

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
  std::size_t sum = 0;

  for (auto [k, info] : seissol::common::enumerate(neighborInfo)) {
    // Need to sort the neighborInfo vectors once
    std::sort(info.second.begin(), info.second.end(), [&](auto a, auto b) {
      return meshTopology.faces()[a].gid() < meshTopology.faces()[b].gid();
    });

    t.insert(info.second.begin(), info.second.end());
    sum += info.second.size();

    // Create MPI neighbor list
    addMPINeighor(meshTopology, info.first, info.second);

    copySide[k].resize(info.second.size());
    ghostSide[k].resize(info.second.size());
    copyFirstVertex[k].resize(info.second.size());
    ghostFirstVertex[k].resize(info.second.size());

    MPI_Irecv(ghostSide[k].data(),
              info.second.size(),
              MPI_CHAR,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[k]);
    MPI_Irecv(ghostFirstVertex[k].data(),
              info.second.size(),
              MPI_UNSIGNED_LONG,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[neighborInfo.size() + k]);

    // Neighbor side
    for (std::size_t i = 0; i < info.second.size(); i++) {
      // The side of boundary
      std::array<int, 2> cellIds{};
      PUML::Upward::cells(meshTopology, faces[info.second[i]], cellIds.data());
      const auto side = PUML::Downward::faceSide(meshTopology, cells[cellIds[0]], info.second[i]);
      logassert(side >= 0 && static_cast<std::size_t>(side) < Cell::NumFaces);
      copySide[k][i] = side;

      std::array<unsigned int, Cell::NumVertices> topoVertices{};
      PUML::Downward::vertices(meshTopology, cells[cellIds[0]], topoVertices.data());

      // First vertex of the face on the boundary
      const auto firstVertex = topoVertices[FirstFaceVertex[PumlFaceToSeisSol[side]]];
      copyFirstVertex[k][i] = vertices[firstVertex].gid();

      // Set the MPI index
      logassert(m_elements[cellIds[0]].mpiIndices[PumlFaceToSeisSol[side]] == 0);
      m_elements[cellIds[0]].mpiIndices[PumlFaceToSeisSol[side]] = i;
    }

    MPI_Isend(copySide[k].data(),
              info.second.size(),
              MPI_CHAR,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[neighborInfo.size() * 2 + k]);
    MPI_Isend(copyFirstVertex[k].data(),
              info.second.size(),
              MPI_UNSIGNED_LONG,
              info.first,
              0,
              Mpi::mpi.comm(),
              &requests[neighborInfo.size() * 3 + k]);
  }
  logassert(t.size() == sum);

  MPI_Waitall(neighborInfo.size() * 4, requests.data(), MPI_STATUSES_IGNORE);

  for (auto [k, info] : seissol::common::enumerate(neighborInfo)) {
    for (std::size_t i = 0; i < info.second.size(); i++) {
      // Set neighbor side
      std::array<int, 2> cellIds{};
      PUML::Upward::cells(meshTopology, faces[info.second[i]], cellIds.data());
      logassert(cellIds[1] < 0);

      // the linters demanded a double cast here
      const auto side = static_cast<std::size_t>(static_cast<unsigned char>(copySide[k][i]));
      const auto gSide = static_cast<std::size_t>(static_cast<unsigned char>(ghostSide[k][i]));
      m_elements[cellIds[0]].neighborSides[PumlFaceToSeisSol[side]] = PumlFaceToSeisSol[gSide];

      // Set side sideOrientation
      std::array<unsigned long, Cell::NumVertices> nvertices{};
      PUML::Downward::gvertices(meshTopology, cells[cellIds[0]], nvertices.data());

      const auto* localFirstVertex =
          std::find(nvertices.begin(), nvertices.end(), ghostFirstVertex[k][i]);
      logassert(localFirstVertex != nvertices.end());

      m_elements[cellIds[0]].sideOrientations[PumlFaceToSeisSol[side]] =
          FaceVertexToOrientation[PumlFaceToSeisSol[side]][localFirstVertex - nvertices.begin()];
      logassert(m_elements[cellIds[0]].sideOrientations[PumlFaceToSeisSol[side]] >= 0);
    }
  }

  // Set vertices
  m_vertices.resize(verticesGeometry.size());
  for (std::size_t i = 0; i < verticesGeometry.size(); i++) {
    memcpy(m_vertices[i].coords, verticesGeometry[i].coordinate(), Cell::Dim * sizeof(double));

    PUML::Upward::cells(meshGeometry, verticesGeometry[i], m_vertices[i].elements);
  }

  // the neighborSide needs to be _inferred_ here.
  for (auto& [_, neighbor] : m_MPINeighbors) {
    for (auto& element : neighbor.elements) {
      element.neighborSide = m_elements[element.localElement].neighborSides[element.localSide];
    }
  }
}

void PUMLReader::addMPINeighor(const PumlMesh& meshTopology,
                               int rank,
                               const std::vector<unsigned int>& faces) {
  const std::size_t id = m_MPINeighbors.size();
  MPINeighbor& neighbor = m_MPINeighbors[rank];

  neighbor.localID = id;
  neighbor.elements.resize(faces.size());

  for (std::size_t i = 0; i < faces.size(); i++) {
    std::array<int, 2> cellIds;
    PUML::Upward::cells(meshTopology, meshTopology.faces()[faces[i]], cellIds.data());

    neighbor.elements[i].localElement = cellIds[0];

    neighbor.elements[i].neighborElement = i;

    std::array<unsigned int, Cell::NumFaces> sides;
    PUML::Downward::faces(meshTopology, meshTopology.cells()[cellIds[0]], sides.data());
    neighbor.elements[i].localSide = [&]() {
      for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
        if (sides[PumlFaceToSeisSol[f]] == faces[i]) {
          return f;
        }
      }
      throw;
    }();
  }
}

bool PUMLReader::inlineTimestepCompute() const { return true; }

bool PUMLReader::inlineClusterCompute() const { return true; }

} // namespace seissol::geometry

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
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Clustering/Clustering.h"
#include "Initializer/Clustering/VertexWeights/VertexWeightModel.h"
#include "Initializer/FaceMap.h"
#include "Initializer/Parameters/MeshParameters.h"
#include "Monitoring/Instrumentation.h"
#include "PartitioningLib.h"

#include <Eigen/Core>
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

namespace seissol::geometry {

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

/**
 * Decodes the boundary condition tag into a string representation.
 */
inline std::string bcToString(uint32_t id, const FaceMap& faceMap) {
  const auto type = faceMap.at(id);
  if (type == FaceType::Regular) {
    return std::string("regular");
  } else if (type == FaceType::FreeSurface) {
    return std::string("free surface");
  } else if (type == FaceType::FreeSurfaceGravity) {
    return std::string("free surface with gravity");
  } else if (type == FaceType::Dirichlet) {
    return std::string("dirichlet");
  } else if (type == FaceType::Outflow) {
    return std::string("outflow");
  } else if (type == FaceType::Analytical) {
    return std::string("analytical");
  } else if (type == FaceType::DynamicRupture) {
    std::stringstream s;
    s << "dynamic rupture (face tag " << id << ")";
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
                                uint8_t side,
                                uint32_t sideBC,
                                uint64_t cellIdAsInFile,
                                const FaceMap& faceMap) {
  // all of these will only issue warnings here -- the "logError()" is supposed to come later, after
  // all warning have been logged

  const auto faceType = faceMap.at(sideBC);

  if (faceType.has_value()) {

    // if a face is an internal face, it has to have a neighbor on either this rank or somewhere
    // else:
    if (getBCType(faceType.value()) == BCType::Internal) {
      if (cellNeighbors[side] < 0 && !face.isShared()) {
        logWarning(true) << "Element" << cellIdAsInFile << ", side" << side << " has a"
                         << bcToString(sideBC, faceMap)
                         << "boundary condition, but the neighboring element doesn't exist";
        return false;
      }
    }
    // external boundaries must not have neighboring elements:
    else {
      if (cellNeighbors[side] >= 0 || face.isShared()) {
        logWarning(true) << "Element" << cellIdAsInFile << ", side" << side << " has a"
                         << bcToString(sideBC, faceMap)
                         << "boundary condition, but a neighboring element exists";
        return false;
      }
    }
  } else {
    // ignore unknown boundary conditions and warn
    logWarning(true) << "Element" << cellIdAsInFile << ", side" << side
                     << " has a boundary condition (" << sideBC
                     << ") which is not understood by this version of SeisSol";
    return false;
  }
  return true;
}

// helper arrays

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

// PUML face p is opposite the local vertex PumlFaceMissingVertex[p]; SeisSol face s is opposite
// the local vertex Cell::NumVertices - 1 - s
const std::array<std::size_t, 4> PumlFaceMissingVertex = {3, 2, 0, 1};

/**
 * Signed volume measure of a tetrahedron, with the same convention as
 * MeshReader::verifyMeshOrientation, which requires a negative value.
 */
double orientationDeterminant(const std::array<const double*, Cell::NumVertices>& coords) {
  Eigen::Matrix<double, 4, 4> mat;
  mat << coords[0][0], coords[0][1], coords[0][2], 1, coords[1][0], coords[1][1], coords[1][2], 1,
      coords[2][0], coords[2][1], coords[2][2], 1, coords[3][0], coords[3][1], coords[3][2], 1;
  return mat.determinant();
}

/**
 * The canonical local vertex numbering of a cell: sort the vertices by topological id, then
 * repair the orientation with the transposition (2 3) where the sorted order would come out
 * positively oriented. Returns the map from the new local slot to the old one.
 *
 * Sorting alone already forces the face orientation index to zero on every interior face.
 * FirstFaceVertex is the smallest local index on each face, hence under a sorted numbering it is
 * the vertex with the smallest topological id on that face -- and both sides of a face agree on
 * which vertex that is. The transposition (2 3) leaves local slots 0 and 1 alone and therefore
 * keeps that property intact, which no other odd permutation does.
 *
 * The sort key has to be the topological id, since that is what identifies the two sides of a
 * face, while the orientation has to be judged on the geometry. For periodic meshes the two
 * differ.
 */
std::array<std::size_t, Cell::NumVertices>
    canonicalVertexOrder(const std::array<unsigned int, Cell::NumVertices>& topoVertices,
                         const std::array<const double*, Cell::NumVertices>& coords) {
  std::array<std::size_t, Cell::NumVertices> order{};
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](auto a, auto b) {
    return topoVertices[a] < topoVertices[b];
  });

  std::array<const double*, Cell::NumVertices> sorted{};
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    sorted[k] = coords[order[k]];
  }
  if (orientationDeterminant(sorted) >= 0) {
    std::swap(order[2], order[3]);
  }
  return order;
}

/**
 * PUML face index -> SeisSol face index for a given canonical vertex order. Both index schemes
 * cover the same four faces; only the ordering within a cell differs.
 */
std::array<std::uint8_t, Cell::NumFaces>
    canonicalFaceMap(const std::array<std::size_t, Cell::NumVertices>& order) {
  std::array<std::size_t, Cell::NumVertices> inverse{};
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    inverse[order[k]] = k;
  }
  std::array<std::uint8_t, Cell::NumFaces> map{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    map[f] = static_cast<std::uint8_t>(Cell::NumVertices - 1 - inverse[PumlFaceMissingVertex[f]]);
  }
  return map;
}

/**
 * Recover the vertex order from a face map. PumlFaceMissingVertex is a bijection, so the face map
 * already determines the permutation and only one of the two needs to be kept per cell.
 */
std::array<std::size_t, Cell::NumVertices>
    vertexOrderFromFaceMap(const std::array<std::uint8_t, Cell::NumFaces>& map) {
  std::array<std::size_t, Cell::NumVertices> order{};
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    order[Cell::NumVertices - 1 - map[f]] = PumlFaceMissingVertex[f];
  }
  return order;
}

template <typename T>
void applyVertexOrder(std::array<T, Cell::NumVertices>& values,
                      const std::array<std::size_t, Cell::NumVertices>& order) {
  const auto original = values;
  for (std::size_t k = 0; k < Cell::NumVertices; ++k) {
    values[k] = original[order[k]];
  }
}
} // namespace

PUMLReader::PUMLReader(const std::string& meshFile,
                       const std::string& partitioningLib,
                       const seissol::FaceMap& faceMap,
                       seissol::initializer::parameters::BoundaryFormat boundaryFormat,
                       seissol::initializer::parameters::TopologyFormat topologyFormat,
                       initializer::Clustering* clustering,
                       initializer::VertexWeightModel* weightModel,
                       double tpwgt) {
  // we need up to two meshes, potentially:
  // one mesh for the geometry
  // one mesh for the topology
  // they will only differ if we have periodic boundary conditions

  PumlMesh meshTopologyExtra;
  PumlMesh meshGeometry;
  meshTopologyExtra.setComm(seissol::Mpi::mpi.comm());
  meshGeometry.setComm(seissol::Mpi::mpi.comm());

  logInfo() << "Read (geometric) connectivity data.";
  read(meshGeometry, meshFile, false, boundaryFormat);

  // Note: we need to call generatePUML in order to create the dual graph of the mesh
  // Note 2: we also need it for vertex identification
  logInfo() << "Generate (geometric) mesh.";
  meshGeometry.generateMesh();

  if (topologyFormat != initializer::parameters::TopologyFormat::Geometric) {
    // we have a topology mesh; separate from the physical mesh

    const bool readTopology =
        topologyFormat == initializer::parameters::TopologyFormat::IdentifyFace;

    logInfo() << "Read (topologic) connectivity data.";
    read(meshTopologyExtra, meshFile, readTopology, boundaryFormat);

    int id = -1;
    if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyVertex) {
      logInfo() << "Read topologic identification data.";
      id = meshTopologyExtra.addData<unsigned long>(
          (std::string(meshFile) + ":/identify").c_str(), PUML::VERTEX, {});
    }

    // generate the topology mesh for the dual graph
    logInfo() << "Generate (topologic) mesh.";
    meshTopologyExtra.generateMesh();

    if (topologyFormat == initializer::parameters::TopologyFormat::IdentifyVertex) {
      // re-identify vertices; then re-distribute
      logInfo() << "Identify topologic vertex data.";
      meshTopologyExtra.identify(id);

      logInfo() << "Generate (topologic) mesh again.";
      meshTopologyExtra.generateMesh();
    }
  }

  auto& meshTopology = topologyFormat == initializer::parameters::TopologyFormat::Geometric
                           ? meshGeometry
                           : meshTopologyExtra;

  // The clustering needs the meshes, which only exist here -- hence the orchestrator is passed
  // in and run rather than its result.
  const initializer::ClusteringResult* clusteringResult = nullptr;
  if (clustering != nullptr) {
    logInfo() << "Compute clustering.";
    clusteringResult = &clustering->compute(meshTopology, meshGeometry);
  }

  logInfo() << "Partition the mesh.";
  partition(meshTopology, meshGeometry, clusteringResult, weightModel, tpwgt, partitioningLib);

  logInfo() << "Generate the correctly-distributed meshes.";
  generatePUML(meshTopology, meshGeometry);

  logInfo() << "Set up mesh data structures.";
  getMesh(meshTopology, meshGeometry, faceMap, boundaryFormat);
}

void PUMLReader::read(PumlMesh& meshTopology,
                      const std::string& file,
                      bool topology,
                      seissol::initializer::parameters::BoundaryFormat boundaryFormat) {
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
                           const initializer::ClusteringResult* clustering,
                           initializer::VertexWeightModel* weightModel,
                           double tpwgt,
                           const std::string& partitioningLib) {
  SCOREP_USER_REGION("PUMLReader_partition", SCOREP_USER_REGION_TYPE_FUNCTION);

  auto partType = toPartitionerType(std::string_view(partitioningLib));
  logInfo() << "Using the" << toStringView(partType) << "partition library and strategy.";
  if (partType == PUML::PartitionerType::None && Mpi::mpi.size() > 1) {
    logWarning()
        << partitioningLib
        << "not found. The performance of this run will probably be much lower than it could be.";
  }
  auto partitioner = PUML::TETPartition::getPartitioner(partType);
  if (partitioner == nullptr) {
    logError() << "Unrecognized partition library: " << partitioningLib;
  }
  auto graph = PUML::TETPartitionGraph(meshTopology);
  weightModel->build(*clustering);
  graph.setVertexWeights(weightModel->vertexWeights(), weightModel->nWeightsPerVertex());

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
  target.setImbalance(weightModel->imbalances()[0] - 1.0);

  auto newPartition = partitioner->partition(graph, target);

  // Written as std::size_t and read back as std::size_t below -- the two have to stay in step,
  // since the read is a reinterpret_cast that would silently misparse a mismatched width.
  meshGeometry.addDataArray(clustering->clusterIds.data(), PUML::CELL, {});
  meshGeometry.addDataArray(clustering->timesteps.cellTimeStepWidths.data(), PUML::CELL, {});

  meshGeometry.partition(newPartition.data());
  if (&meshTopology != &meshGeometry) {
    meshTopology.partition(newPartition.data());
  }
}

void PUMLReader::generatePUML(PumlMesh& meshTopology, PumlMesh& meshGeometry) {
  SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

  if (&meshTopology != &meshGeometry) {
    logInfo() << "Generate the correct (topologic) mesh.";
    meshTopology.generateMesh();
  }
  logInfo() << "Generate the correct (geometric) mesh.";
  meshGeometry.generateMesh();
}

void PUMLReader::getMesh(const PumlMesh& meshTopology,
                         const PumlMesh& meshGeometry,
                         const FaceMap& faceMap,
                         seissol::initializer::parameters::BoundaryFormat boundaryFormat) {
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
  const auto* clusterIds = reinterpret_cast<const std::size_t*>(meshGeometry.cellData(3));
  const auto* timestep = reinterpret_cast<const double*>(meshGeometry.cellData(4));

  std::unordered_map<int, std::vector<unsigned int>> neighborInfo; // List of shared local face ids

  bool isMeshCorrect = true;

  // Canonical local vertex numbering. Computed for every cell up front, because the neighbor
  // lookup below needs the numbering of cells that come later in the loop.
  std::vector<std::array<std::uint8_t, Cell::NumFaces>> pumlFaceMaps(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    std::array<unsigned int, Cell::NumVertices> topoVertices{};
    PUML::Downward::vertices(meshTopology, cells[i], topoVertices.data());

    std::array<unsigned int, Cell::NumVertices> geomVertices{};
    PUML::Downward::vertices(meshGeometry, cellsGeometry[i], geomVertices.data());

    std::array<const double*, Cell::NumVertices> coords{};
    for (std::size_t k = 0; k < Cell::NumVertices; k++) {
      coords[k] = verticesGeometry[geomVertices[k]].coordinate();
    }

    pumlFaceMaps[i] = canonicalFaceMap(canonicalVertexOrder(topoVertices, coords));
  }

  // Compute everything local
  elements_.resize(cells.size());
  for (std::size_t i = 0; i < cells.size(); i++) {
    elements_[i].globalId = cellIdsAsInFile[i];
    elements_[i].localId = i;
    elements_[i].clusterId = clusterIds[i];
    elements_[i].timestep = timestep[i];

    const auto& pumlToSeisSol = pumlFaceMaps[i];
    const auto vertexOrder = vertexOrderFromFaceMap(pumlToSeisSol);

    // Vertices
    std::array<unsigned int, Cell::NumVertices> geomVertices{};
    PUML::Downward::vertices(meshGeometry, cellsGeometry[i], geomVertices.data());
    applyVertexOrder(geomVertices, vertexOrder);
    for (std::size_t k = 0; k < Cell::NumVertices; k++) {
      elements_[i].vertices[k] = geomVertices[k];
    }

    std::array<unsigned int, Cell::NumVertices> topoVertices{};
    PUML::Downward::vertices(meshTopology, cells[i], topoVertices.data());
    applyVertexOrder(topoVertices, vertexOrder);

    // Neighbor information
    std::array<unsigned int, Cell::NumFaces> faceids{};
    PUML::Downward::faces(meshTopology, cells[i], faceids.data());
    std::array<int, Cell::NumFaces> neighbors{};
    PUML::Neighbor::face(meshTopology, i, neighbors.data());

    for (std::size_t j = 0; j < Cell::NumFaces; j++) {
      const auto faceTag = decodeBoundary(boundaryCond, i, j, boundaryFormat);
      const bool isLocallyCorrect = checkMeshCorrectnessLocally<PumlTopology>(
          faces[faceids[j]], neighbors, j, faceTag, cellIdsAsInFile[i], faceMap);
      isMeshCorrect &= isLocallyCorrect;
      const auto side = pumlToSeisSol[j];

      if (neighbors[j] < 0) {
        elements_[i].neighbors[side] = cellsGeometry.size();

        if (!faces[faceids[j]].isShared()) {
          // Boundary sides
          elements_[i].neighborRanks[side] = rank;
        } else {
          // MPI Boundary
          neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

          elements_[i].neighborRanks[side] = faces[faceids[j]].shared()[0];
        }
      } else {
        logassert(neighbors[j] >= 0 && static_cast<std::size_t>(neighbors[j]) < cells.size());

        elements_[i].neighbors[side] = neighbors[j];

        std::array<int, Cell::NumFaces> nfaces{};
        PUML::Neighbor::face(meshTopology, neighbors[j], nfaces.data());
        const auto* back = std::find(nfaces.begin(), nfaces.end(), i);
        logassert(back != nfaces.end());

        const auto& neighborPumlToSeisSol = pumlFaceMaps[neighbors[j]];
        elements_[i].neighborSides[side] = neighborPumlToSeisSol[back - nfaces.begin()];

        const auto firstVertex = topoVertices[FirstFaceVertex[side]];

        std::array<unsigned int, Cell::NumVertices> nvertices{};
        PUML::Downward::vertices(meshTopology, cells[neighbors[j]], nvertices.data());
        applyVertexOrder(nvertices, vertexOrderFromFaceMap(neighborPumlToSeisSol));
        const auto* neighborFirstVertex =
            std::find(nvertices.begin(), nvertices.end(), firstVertex);
        logassert(neighborFirstVertex != nvertices.end());

        elements_[i].sideOrientations[side] =
            FaceVertexToOrientation[elements_[i].neighborSides[side]]
                                   [neighborFirstVertex - nvertices.begin()];
        logassert(elements_[i].sideOrientations[side] == 0);

        elements_[i].neighborRanks[side] = rank;
      }

      const auto bcCurrentFace = faceMap.at(faceTag);

      elements_[i].boundaries[side] = bcCurrentFace.value_or(FaceType::Regular);
      elements_[i].faultTags[side] = faceTag;
      elements_[i].mpiIndices[side] = 0;
    }

    elements_[i].group = group[i];
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
    addMPINeighor(meshTopology, info.first, info.second, pumlFaceMaps);

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
      const auto pumlSide =
          PUML::Downward::faceSide(meshTopology, cells[cellIds[0]], info.second[i]);
      logassert(pumlSide >= 0 && static_cast<std::size_t>(pumlSide) < Cell::NumFaces);

      // Send the SeisSol side under the canonical numbering: the receiving rank cannot translate
      // a PUML side, because the numbering is per cell now.
      const auto vertexOrder = vertexOrderFromFaceMap(pumlFaceMaps[cellIds[0]]);
      const auto side = pumlFaceMaps[cellIds[0]][pumlSide];
      copySide[k][i] = static_cast<char>(side);

      std::array<unsigned int, Cell::NumVertices> topoVertices{};
      PUML::Downward::vertices(meshTopology, cells[cellIds[0]], topoVertices.data());
      applyVertexOrder(topoVertices, vertexOrder);

      // First vertex of the face on the boundary
      const auto firstVertex = topoVertices[FirstFaceVertex[side]];
      copyFirstVertex[k][i] = vertices[firstVertex].gid();

      // Set the MPI index
      logassert(elements_[cellIds[0]].mpiIndices[side] == 0);
      elements_[cellIds[0]].mpiIndices[side] = i;
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

      // the linters demanded a double cast here; both sides are already SeisSol sides under the
      // canonical numbering, so no translation is applied to either of them
      const auto side = static_cast<std::size_t>(static_cast<unsigned char>(copySide[k][i]));
      const auto gSide = static_cast<std::size_t>(static_cast<unsigned char>(ghostSide[k][i]));
      elements_[cellIds[0]].neighborSides[side] = gSide;

      // Set side sideOrientation
      std::array<unsigned long, Cell::NumVertices> nvertices{};
      PUML::Downward::gvertices(meshTopology, cells[cellIds[0]], nvertices.data());
      applyVertexOrder(nvertices, vertexOrderFromFaceMap(pumlFaceMaps[cellIds[0]]));

      const auto* localFirstVertex =
          std::find(nvertices.begin(), nvertices.end(), ghostFirstVertex[k][i]);
      logassert(localFirstVertex != nvertices.end());

      elements_[cellIds[0]].sideOrientations[side] =
          FaceVertexToOrientation[side][localFirstVertex - nvertices.begin()];
      logassert(elements_[cellIds[0]].sideOrientations[side] == 0);
    }
  }

  // Set vertices
  vertices_.resize(verticesGeometry.size());
  for (std::size_t i = 0; i < verticesGeometry.size(); i++) {
    memcpy(vertices_[i].coords, verticesGeometry[i].coordinate(), Cell::Dim * sizeof(double));

    PUML::Upward::cells(meshGeometry, verticesGeometry[i], vertices_[i].elements);
  }

  // the neighborSide needs to be _inferred_ here.
  for (auto& [_, neighbor] : mpiNeighbors_) {
    for (auto& element : neighbor.elements) {
      element.neighborSide = elements_[element.localElement].neighborSides[element.localSide];
    }
  }
}

void PUMLReader::addMPINeighor(
    const PumlMesh& meshTopology,
    int rank,
    const std::vector<unsigned int>& faces,
    const std::vector<std::array<std::uint8_t, Cell::NumFaces>>& pumlFaceMaps) {
  const std::size_t id = mpiNeighbors_.size();
  MPINeighbor& neighbor = mpiNeighbors_[rank];

  neighbor.localID = id;
  neighbor.elements.resize(faces.size());

  for (std::size_t i = 0; i < faces.size(); i++) {
    std::array<int, 2> cellIds{};
    PUML::Upward::cells(meshTopology, meshTopology.faces()[faces[i]], cellIds.data());

    neighbor.elements[i].localElement = cellIds[0];

    neighbor.elements[i].neighborElement = i;

    std::array<unsigned int, Cell::NumFaces> sides{};
    PUML::Downward::faces(meshTopology, meshTopology.cells()[cellIds[0]], sides.data());
    neighbor.elements[i].localSide = [&]() -> std::size_t {
      for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
        if (sides[f] == faces[i]) {
          return pumlFaceMaps[cellIds[0]][f];
        }
      }
      throw;
    }();
  }
}

bool PUMLReader::inlineTimestepCompute() const { return true; }

bool PUMLReader::inlineClusterCompute() const { return true; }

} // namespace seissol::geometry

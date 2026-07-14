// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MeshReader.h"

#include "Common/Constants.h"
#include "Common/Iterator.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/ParameterDB.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/TimeStepping/GlobalTimestep.h"
#include "MeshDefinition.h"
#include "MeshTools.h"
#include "Parallel/MPI.h"
#include "SeisSol.h"

#include <Eigen/Core>
#include <PUML/TypeInference.h>
#include <algorithm>
#include <cstddef>
#include <map>
#include <mpi.h>
#include <unordered_map>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::geometry {

MeshReader::MeshReader() = default;

MeshReader::~MeshReader() = default;

const std::vector<Element>& MeshReader::getElements() const { return elements_; }

const std::vector<Vertex>& MeshReader::getVertices() const { return vertices_; }

const std::map<int, MPINeighbor>& MeshReader::getMPINeighbors() const { return mpiNeighbors_; }

const std::map<int, std::vector<MPINeighborElement>>& MeshReader::getMPIFaultNeighbors() const {
  return mpiFaultNeighbors_;
}

const std::unordered_map<int, std::vector<GhostElementMetadata>>&
    MeshReader::getGhostlayerMetadata() const {
  return ghostlayerMetadata_;
}

const std::vector<Fault>& MeshReader::getFault() const { return fault_; }

bool MeshReader::hasFault() const { return !fault_.empty(); }

bool MeshReader::hasPlusFault() const { return hasPlusFault_; }

const std::vector<LinearGhostCell>& MeshReader::linearGhostlayer() const {
  return linearGhostlayer_;
}

const std::map<std::pair<int, std::size_t>, std::size_t>& MeshReader::toLinearGhostlayer() const {
  return toLinearGhostlayer_;
}

void MeshReader::displaceMesh(const Eigen::Vector3d& displacement) {
  for (std::size_t vertexNo = 0; vertexNo < vertices_.size(); ++vertexNo) {
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      vertices_[vertexNo].coords[i] += displacement[i];
    }
  }
}

// TODO: Test proper scaling
//  scalingMatrix is stored column-major, i.e.
//  scalingMatrix_ij = scalingMatrix[j][i]
void MeshReader::scaleMesh(const Eigen::Matrix3d& scalingMatrix) {
  for (std::size_t vertexNo = 0; vertexNo < vertices_.size(); ++vertexNo) {
    Eigen::Vector3d point;
    point << vertices_[vertexNo].coords[0], vertices_[vertexNo].coords[1],
        vertices_[vertexNo].coords[2];
    const auto result = scalingMatrix * point;
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      vertices_[vertexNo].coords[i] = result[i];
    }
  }
}

void MeshReader::disableDR() {
  for (auto& elem : elements_) {
    for (std::size_t j = 0; j < Cell::NumFaces; ++j) {
      if (elem.boundaries[j] == FaceType::DynamicRupture) {
        elem.boundaries[j] = FaceType::Regular;
      }
    }
  }
}

/**
 * Reconstruct the fault information from the boundary conditions
 */
void MeshReader::extractFaultInformation(
    const VrtxCoords& refPoint, seissol::initializer::parameters::RefPointMethod refPointMethod) {
  for (auto& cell : elements_) {

    for (std::size_t j = 0; j < Cell::NumFaces; ++j) {
      // Set default mpi fault indices
      cell.mpiFaultIndices[j] = -1;

      if (cell.boundaries[j] != FaceType::DynamicRupture) {
        continue;
      }

      // DR boundary

      if (cell.neighborRanks[j] == seissol::Mpi::mpi.rank()) {
        // Completely local DR boundary

        if (cell.neighbors[j] < cell.localId) {
          // This was already handled by the other side
          continue;
        }
      } else {
        // Handle boundary faces

        // FIXME we use the MPI number here for the neighbor element id
        // It is not very nice but should generate the correct ordering.
        const MPINeighborElement neighbor = {
            cell.localId, static_cast<SideId>(j), cell.mpiIndices[j], cell.neighborSides[j]};
        mpiFaultNeighbors_[cell.neighborRanks[j]].push_back(neighbor);
      }

      Fault f{};

      // Detect +/- side

      // Compute normal of the DR face
      // Boundary side vector pointing in chi- and tau-direction
      VrtxCoords chiVec;
      VrtxCoords tauVec;
      MeshTools::sub(vertices_[cell.vertices[MeshTools::FACE2NODES[j][1]]].coords,
                     vertices_[cell.vertices[MeshTools::FACE2NODES[j][0]]].coords,
                     chiVec);
      MeshTools::sub(vertices_[cell.vertices[MeshTools::FACE2NODES[j][2]]].coords,
                     vertices_[cell.vertices[MeshTools::FACE2NODES[j][0]]].coords,
                     tauVec);
      MeshTools::cross(chiVec, tauVec, f.normal);

      // Normalize normal
      MeshTools::mul(f.normal, 1.0 / MeshTools::norm(f.normal), f.normal);

      // Check whether the tetrahedron and the reference point are on the same side of the face
      VrtxCoords tmp1;
      VrtxCoords tmp2;
      MeshTools::sub(refPoint, vertices_[cell.vertices[MeshTools::FACE2NODES[j][0]]].coords, tmp1);
      MeshTools::sub(vertices_[cell.vertices[MeshTools::FACE2MISSINGNODE[j]]].coords,
                     vertices_[cell.vertices[MeshTools::FACE2NODES[j][0]]].coords,
                     tmp2);
      bool isPlus = false;
      if (refPointMethod == seissol::initializer::parameters::RefPointMethod::Point) {
        isPlus = MeshTools::dot(tmp1, f.normal) * MeshTools::dot(tmp2, f.normal) > 0;
      } else {
        isPlus = MeshTools::dot(refPoint, f.normal) > 0;
      }
      // Fix normal direction and get correct chiVec
      if (!isPlus) {
        // In case of a minus side, compute chi using node 0 and 1 from the plus side
        MeshTools::sub(
            vertices_[cell.vertices
                          [MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                        [(3 + 1 - cell.sideOrientations[j]) % 3]]]]
                .coords,
            vertices_[cell.vertices
                          [MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                        [(3 + 0 - cell.sideOrientations[j]) % 3]]]]
                .coords,
            chiVec);

        MeshTools::sub(
            vertices_[cell.vertices
                          [MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                        [(3 + 2 - cell.sideOrientations[j]) % 3]]]]
                .coords,
            vertices_[cell.vertices
                          [MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                        [(3 + 0 - cell.sideOrientations[j]) % 3]]]]
                .coords,
            tauVec);

        MeshTools::cross(chiVec, tauVec, f.normal);
        MeshTools::mul(f.normal, 1.0 / MeshTools::norm(f.normal), f.normal);
      }

      // Compute vector inside the triangle's plane for the rotation matrix
      MeshTools::mul(chiVec, 1.0 / MeshTools::norm(chiVec), f.tangent1);
      // Compute second vector in the plane, orthogonal to the normal and tangent 1 vectors
      MeshTools::cross(f.normal, f.tangent1, f.tangent2);

      auto remoteNeighbor = cell.neighbors[j] == static_cast<int>(elements_.size());

      // Index of the element on the other side
      const int neighborIndex = remoteNeighbor ? -1 : cell.neighbors[j];

      const GlobalElemId neighborGlobalId =
          remoteNeighbor ? ghostlayerMetadata_[cell.neighborRanks[j]][cell.mpiIndices[j]].globalId
                         : elements_[cell.neighbors[j]].globalId;

      if (isPlus) {
        f.globalId = cell.globalId;
        f.element = cell.localId;
        f.side = j;
        f.neighborGlobalId = neighborGlobalId;
        f.neighborElement = neighborIndex;
        f.neighborSide = cell.neighborSides[j];
        f.tag = cell.faultTags[j];
      } else {
        f.globalId = neighborGlobalId;
        f.element = neighborIndex;
        f.side = cell.neighborSides[j];
        f.neighborGlobalId = cell.globalId;
        f.neighborElement = cell.localId;
        f.neighborSide = j;
        f.tag = cell.faultTags[j];
      }

      fault_.push_back(f);

      // Check if we have a plus fault side
      if (isPlus || neighborIndex >= 0) {
        hasPlusFault_ = true;
      }
    }
  }

  // Sort fault neighbor lists and update MPI fault indices
  for (auto& neighbor : mpiFaultNeighbors_) {

    if (neighbor.first > seissol::Mpi::mpi.rank()) {
      std::sort(neighbor.second.begin(),
                neighbor.second.end(),
                [](const MPINeighborElement& elem1, const MPINeighborElement& elem2) {
                  return (elem1.localElement < elem2.localElement) ||
                         (elem1.localElement == elem2.localElement &&
                          elem1.localSide < elem2.localSide);
                });
    } else {
      std::sort(neighbor.second.begin(),
                neighbor.second.end(),
                [](const MPINeighborElement& elem1, const MPINeighborElement& elem2) {
                  return (elem1.neighborElement < elem2.neighborElement) ||
                         (elem1.neighborElement == elem2.neighborElement &&
                          elem1.neighborSide < elem2.neighborSide);
                });
    }

    // Set the MPI fault number of all elements
    for (std::size_t j = 0; j < neighbor.second.size(); ++j) {
      elements_[neighbor.second[j].localElement].mpiFaultIndices[neighbor.second[j].localSide] = j;
    }
  }
}

void MeshReader::exchangeGhostlayerMetadata() {
  std::unordered_map<int, std::vector<GhostElementMetadata>> sendData;
  std::unordered_map<int, std::vector<GhostElementMetadata>> recvData;

  constexpr int Tag = 10;
  MPI_Comm comm = seissol::Mpi::mpi.comm();

  std::vector<MPI_Request> requests(mpiNeighbors_.size() * 2);

  // TODO(David): Once a generic MPI type inference module is ready, replace this part here ...
  // Maybe.
  MPI_Datatype ghostElementType = MPI_DATATYPE_NULL;
  MPI_Datatype ghostElementTypePre = MPI_DATATYPE_NULL;

  // assume that all vertices are stored contiguously
  const int datatypeCount = 6;
  const std::vector<int> datatypeBlocklen{Cell::NumVertices * Cell::Dim, 1, 1, 1, 1, 1};
  const std::vector<MPI_Aint> datatypeDisplacement{offsetof(GhostElementMetadata, vertices),
                                                   offsetof(GhostElementMetadata, group),
                                                   offsetof(GhostElementMetadata, localId),
                                                   offsetof(GhostElementMetadata, globalId),
                                                   offsetof(GhostElementMetadata, clusterId),
                                                   offsetof(GhostElementMetadata, timestep)};
  const std::vector<MPI_Datatype> datatypeDatatype{MPI_DOUBLE,
                                                   MPI_INT,
                                                   PUML::MPITypeInfer<LocalElemId>::type(),
                                                   PUML::MPITypeInfer<GlobalElemId>::type(),
                                                   MPI_INT,
                                                   MPI_DOUBLE};

  MPI_Type_create_struct(datatypeCount,
                         datatypeBlocklen.data(),
                         datatypeDisplacement.data(),
                         datatypeDatatype.data(),
                         &ghostElementTypePre);
  MPI_Type_create_resized(ghostElementTypePre, 0, sizeof(GhostElementMetadata), &ghostElementType);
  MPI_Type_commit(&ghostElementType);

  size_t counter = 0;
  for (auto it = mpiNeighbors_.begin(); it != mpiNeighbors_.end(); ++it, counter += 2) {
    const auto targetRank = it->first;
    const auto count = it->second.elements.size();

    recvData[targetRank].resize(count);

    sendData[targetRank].resize(count);
    for (size_t j = 0; j < count; ++j) {
      const auto elementIdx = it->second.elements[j].localElement;
      const auto& element = elements_.at(elementIdx);
      auto& ghost = sendData[targetRank][j];

      for (size_t v = 0; v < Cell::NumVertices; ++v) {
        const auto& vertex = vertices_[element.vertices[v]];
        for (std::size_t d = 0; d < Cell::Dim; ++d) {
          ghost.vertices[v][d] = vertex.coords[d];
        }
      }
      ghost.group = element.group;
      ghost.localId = element.localId;
      ghost.globalId = element.globalId;
      ghost.clusterId = element.clusterId;
      ghost.timestep = element.timestep;
    }

    // TODO(David): evaluate, if MPI_Ssend (instead of just MPI_Send) makes sense here?
    MPI_Irecv(recvData[targetRank].data(),
              count,
              ghostElementType,
              targetRank,
              Tag,
              comm,
              &requests[counter]);
    MPI_Isend(sendData[targetRank].data(),
              count,
              ghostElementType,
              targetRank,
              Tag,
              comm,
              &requests[counter + 1]);
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

  ghostlayerMetadata_ = std::move(recvData);

  MPI_Type_free(&ghostElementType);
}

void MeshReader::linearizeGhostlayer() {
  linearGhostlayer_.clear();
  toLinearGhostlayer_.clear();

  // basic assumption: each cell appears only on exactly one rank
  for (const auto& [rank, cells] : ghostlayerMetadata_) {
    // map for deterministic ordering (for now)
    std::map<std::size_t, std::vector<std::size_t>> ordering;
    for (std::size_t i = 0; i < cells.size(); ++i) {
      ordering[cells[i].globalId].emplace_back(i);
    }
    for (const auto& [_, ids] : ordering) {
      for (const auto& index : ids) {
        toLinearGhostlayer_[{rank, index}] = linearGhostlayer_.size();
      }
      linearGhostlayer_.push_back(LinearGhostCell{ids, rank});
    }
  }
}

void MeshReader::computeTimestepIfNecessary(const seissol::SeisSol& seissolInstance) {
  if (!inlineTimestepCompute()) {
    const auto ctvarray = seissol::initializer::CellToVertexArray::fromMeshReader(*this);
    const auto timesteps =
        seissol::initializer::computeTimesteps(ctvarray, seissolInstance.parameters());
    for (auto [cell, timestep] : seissol::common::zip(elements_, timesteps.cellTimeStepWidths)) {
      cell.timestep = timestep;

      // enforce GTS in the case here
      cell.clusterId = 0;
    }
  }
}

void MeshReader::verifyMeshOrientation() {
  // for now, only check the tetrahedron orientation here

  bool correct = true;

#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < elements_.size(); ++i) {
    const auto& element = elements_[i];

    // check orientation
    Eigen::Matrix<double, 4, 4> mat;
    const auto& v1 = vertices_[element.vertices[0]].coords;
    const auto& v2 = vertices_[element.vertices[1]].coords;
    const auto& v3 = vertices_[element.vertices[2]].coords;
    const auto& v4 = vertices_[element.vertices[3]].coords;

    mat << v1[0], v1[1], v1[2], 1, v2[0], v2[1], v2[2], 1, v3[0], v3[1], v3[2], 1, v4[0], v4[1],
        v4[2], 1;

    if (mat.determinant() >= 0) {
      logWarning()
          << "The cell" << i
          << "has an incorrect orientation. Thus, SeisSol would produce incorrect results.";
      correct = false;
    }
  }

  if (!correct) {
    logError() << "There are geometric problems with the given mesh.";
  }
}

} // namespace seissol::geometry

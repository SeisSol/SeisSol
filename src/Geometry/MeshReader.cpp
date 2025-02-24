// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MeshReader.h"

#include "MeshDefinition.h"
#include "MeshTools.h"

#include "PUML/TypeInference.h"
#include "Parallel/MPI.h"
#include <Common/Iterator.h>
#include <Initializer/ParameterDB.h>
#include <Initializer/Parameters/DRParameters.h>
#include <Initializer/TimeStepping/GlobalTimestep.h>
#include <SeisSol.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif

namespace seissol::geometry {

MeshReader::MeshReader(int rank) : mRank(rank) {}

MeshReader::~MeshReader() = default;

const std::vector<Element>& MeshReader::getElements() const { return m_elements; }

const std::vector<Vertex>& MeshReader::getVertices() const { return m_vertices; }

const std::map<int, MPINeighbor>& MeshReader::getMPINeighbors() const { return m_MPINeighbors; }

const std::map<int, std::vector<MPINeighborElement>>& MeshReader::getMPIFaultNeighbors() const {
  return m_MPIFaultNeighbors;
}

const std::unordered_map<int, std::vector<GhostElementMetadata>>&
    MeshReader::getGhostlayerMetadata() const {
  return m_ghostlayerMetadata;
}

const std::vector<Fault>& MeshReader::getFault() const { return m_fault; }

bool MeshReader::hasFault() const { return !m_fault.empty(); }

bool MeshReader::hasPlusFault() const { return m_hasPlusFault; }

void MeshReader::displaceMesh(const Eigen::Vector3d& displacement) {
  for (unsigned vertexNo = 0; vertexNo < m_vertices.size(); ++vertexNo) {
    for (unsigned i = 0; i < 3; ++i) {
      m_vertices[vertexNo].coords[i] += displacement[i];
    }
  }
}

// TODO: Test proper scaling
//  scalingMatrix is stored column-major, i.e.
//  scalingMatrix_ij = scalingMatrix[j][i]
void MeshReader::scaleMesh(const Eigen::Matrix3d& scalingMatrix) {
  for (unsigned vertexNo = 0; vertexNo < m_vertices.size(); ++vertexNo) {
    Eigen::Vector3d point;
    point << m_vertices[vertexNo].coords[0], m_vertices[vertexNo].coords[1],
        m_vertices[vertexNo].coords[2];
    const auto result = scalingMatrix * point;
    for (unsigned i = 0; i < 3; ++i) {
      m_vertices[vertexNo].coords[i] = result[i];
    }
  }
}

/**
 * Reconstruct the fault information from the boundary conditions
 */
void MeshReader::extractFaultInformation(
    const VrtxCoords& refPoint, seissol::initializer::parameters::RefPointMethod refPointMethod) {
  for (auto& i : m_elements) {

    for (int j = 0; j < 4; j++) {
      // Set default mpi fault indices
      i.mpiFaultIndices[j] = -1;

      if (i.boundaries[j] != 3) {
        continue;
      }

      // DR boundary

      if (i.neighborRanks[j] == mRank) {
        // Completely local DR boundary

        if (i.neighbors[j] < i.localId) {
          // This was already handled by the other side
          continue;
        }
      } else {
        // Handle boundary faces

        // FIXME we use the MPI number here for the neighbor element id
        // It is not very nice but should generate the correct ordering.
        const MPINeighborElement neighbor = {i.localId, j, i.mpiIndices[j], i.neighborSides[j]};
        m_MPIFaultNeighbors[i.neighborRanks[j]].push_back(neighbor);
      }

      Fault f{};

      // Detect +/- side
      // Computes the distance between the bary center of the tetrahedron and the face
      // Does not work for all meshes
      //                VrtxCoords elementCenter;
      //                VrtxCoords faceCenter;
      //                MeshTools::center(*i, m_vertices, elementCenter);
      //                MeshTools::center(*i, j, m_vertices, faceCenter);
      //
      //                bool isPlus = (MeshTools::distance(elementCenter, center)
      //                        < MeshTools::distance(faceCenter, center));

      // Compute normal of the DR face
      // Boundary side vector pointing in chi- and tau-direction
      VrtxCoords chiVec;
      VrtxCoords tauVec;
      MeshTools::sub(m_vertices[i.vertices[MeshTools::FACE2NODES[j][1]]].coords,
                     m_vertices[i.vertices[MeshTools::FACE2NODES[j][0]]].coords,
                     chiVec);
      MeshTools::sub(m_vertices[i.vertices[MeshTools::FACE2NODES[j][2]]].coords,
                     m_vertices[i.vertices[MeshTools::FACE2NODES[j][0]]].coords,
                     tauVec);
      MeshTools::cross(chiVec, tauVec, f.normal);

      // Normalize normal
      MeshTools::mul(f.normal, 1.0 / MeshTools::norm(f.normal), f.normal);

      // Check whether the tetrahedron and the reference point are on the same side of the face
      VrtxCoords tmp1;
      VrtxCoords tmp2;
      MeshTools::sub(refPoint, m_vertices[i.vertices[MeshTools::FACE2NODES[j][0]]].coords, tmp1);
      MeshTools::sub(m_vertices[i.vertices[MeshTools::FACE2MISSINGNODE[j]]].coords,
                     m_vertices[i.vertices[MeshTools::FACE2NODES[j][0]]].coords,
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
            m_vertices
                [i.vertices[MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                         [(3 + 1 - i.sideOrientations[j]) % 3]]]]
                    .coords,
            m_vertices
                [i.vertices[MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                         [(3 + 0 - i.sideOrientations[j]) % 3]]]]
                    .coords,
            chiVec);

        MeshTools::sub(
            m_vertices
                [i.vertices[MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                         [(3 + 2 - i.sideOrientations[j]) % 3]]]]
                    .coords,
            m_vertices
                [i.vertices[MeshTools::FACE2NODES[j][MeshTools::NEIGHBORFACENODE2LOCAL
                                                         [(3 + 0 - i.sideOrientations[j]) % 3]]]]
                    .coords,
            tauVec);

        MeshTools::cross(chiVec, tauVec, f.normal);
        MeshTools::mul(f.normal, 1.0 / MeshTools::norm(f.normal), f.normal);
      }

      // Compute vector inside the triangle's plane for the rotation matrix
      MeshTools::mul(chiVec, 1.0 / MeshTools::norm(chiVec), f.tangent1);
      // Compute second vector in the plane, orthogonal to the normal and tangent 1 vectors
      MeshTools::cross(f.normal, f.tangent1, f.tangent2);

      auto remoteNeighbor = i.neighbors[j] == static_cast<int>(m_elements.size());

      // Index of the element on the other side
      const int neighborIndex = remoteNeighbor ? -1 : i.neighbors[j];

      const GlobalElemId neighborGlobalId =
          remoteNeighbor ? m_ghostlayerMetadata[i.neighborRanks[j]][i.mpiIndices[j]].globalId
                         : m_elements[i.neighbors[j]].globalId;

      if (isPlus) {
        f.globalId = i.globalId;
        f.element = i.localId;
        f.side = j;
        f.neighborGlobalId = neighborGlobalId;
        f.neighborElement = neighborIndex;
        f.neighborSide = i.neighborSides[j];
        f.tag = i.faultTags[j];
      } else {
        f.globalId = neighborGlobalId;
        f.element = neighborIndex;
        f.side = i.neighborSides[j];
        f.neighborGlobalId = i.globalId;
        f.neighborElement = i.localId;
        f.neighborSide = j;
        f.tag = i.faultTags[j];
      }

      m_fault.push_back(f);

      // Check if we have a plus fault side
      if (isPlus || neighborIndex >= 0) {
        m_hasPlusFault = true;
      }
    }
  }

  // Sort fault neighbor lists and update MPI fault indices
  for (auto& i : m_MPIFaultNeighbors) {

    if (i.first > mRank) {
      std::sort(i.second.begin(),
                i.second.end(),
                [](const MPINeighborElement& elem1, const MPINeighborElement& elem2) {
                  return (elem1.localElement < elem2.localElement) ||
                         (elem1.localElement == elem2.localElement &&
                          elem1.localSide < elem2.localSide);
                });
    } else {
      std::sort(i.second.begin(),
                i.second.end(),
                [](const MPINeighborElement& elem1, const MPINeighborElement& elem2) {
                  return (elem1.neighborElement < elem2.neighborElement) ||
                         (elem1.neighborElement == elem2.neighborElement &&
                          elem1.neighborSide < elem2.neighborSide);
                });
    }

    // Set the MPI fault number of all elements
    for (int j = 0; j < static_cast<int>(i.second.size()); j++) {
      m_elements[i.second[j].localElement].mpiFaultIndices[i.second[j].localSide] = j;
    }
  }
}

void MeshReader::exchangeGhostlayerMetadata() {
#ifdef USE_MPI
  std::unordered_map<int, std::vector<GhostElementMetadata>> sendData;
  std::unordered_map<int, std::vector<GhostElementMetadata>> recvData;

  constexpr int Tag = 10;
  MPI_Comm comm = seissol::MPI::mpi.comm();

  std::vector<MPI_Request> requests(m_MPINeighbors.size() * 2);

  // TODO(David): Once a generic MPI type inference module is ready, replace this part here ...
  // Maybe.
  MPI_Datatype ghostElementType = MPI_DATATYPE_NULL;

  // assume that all vertices are stored contiguously
  const int datatypeCount = 6;
  const std::vector<int> datatypeBlocklen{12, 1, 1, 1, 1, 1};
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
                         &ghostElementType);
  MPI_Type_commit(&ghostElementType);

  size_t counter = 0;
  for (auto it = m_MPINeighbors.begin(); it != m_MPINeighbors.end(); ++it, counter += 2) {
    const auto targetRank = it->first;
    const auto count = it->second.elements.size();

    recvData[targetRank].resize(count);

    sendData[targetRank].resize(count);
    for (size_t j = 0; j < count; ++j) {
      const auto elementIdx = it->second.elements[j].localElement;
      const auto& element = m_elements.at(elementIdx);
      auto& ghost = sendData[targetRank][j];

      for (size_t v = 0; v < 4; ++v) {
        const auto& vertex = m_vertices[element.vertices[v]];
        ghost.vertices[v][0] = vertex.coords[0];
        ghost.vertices[v][1] = vertex.coords[1];
        ghost.vertices[v][2] = vertex.coords[2];
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

  m_ghostlayerMetadata = std::move(recvData);

  MPI_Type_free(&ghostElementType);
#endif
}

void MeshReader::computeTimestepIfNecessary(const seissol::SeisSol& seissolInstance) {
  if (!inlineTimestepCompute()) {
    const auto ctvarray = seissol::initializer::CellToVertexArray::fromMeshReader(*this);
    const auto timesteps =
        seissol::initializer::computeTimesteps(ctvarray, seissolInstance.getSeisSolParameters());
    for (auto [cell, timestep] : seissol::common::zip(m_elements, timesteps.cellTimeStepWidths)) {
      cell.timestep = timestep;

      // enforce GTS in the case here
      cell.clusterId = 0;
    }
  }
}

} // namespace seissol::geometry

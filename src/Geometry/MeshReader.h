// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_MESHREADER_H_
#define SEISSOL_SRC_GEOMETRY_MESHREADER_H_

#include "MeshDefinition.h"

#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "Initializer/Parameters/DRParameters.h"

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::geometry {

struct GhostElementMetadata {
  double vertices[4][3];
  int group;
  LocalElemId localId;
  GlobalElemId globalId;
  int clusterId;
  double timestep;
};

class MeshReader {
  protected:
  const int mRank;

  std::vector<Element> m_elements;

  std::vector<Vertex> m_vertices;

  /** Convert global element index to local */
  std::map<int, int> m_g2lElements;

  /** Convert global vertex index to local */
  std::map<int, int> m_g2lVertices;

  /** Number of MPI neighbors */
  std::map<int, MPINeighbor> m_MPINeighbors;

  /** Number of MPI fault neighbors */
  std::map<int, std::vector<MPINeighborElement>> m_MPIFaultNeighbors;

  /** Fault information */
  std::vector<Fault> m_fault;

  /** Vertices of MPI Neighbors*/
  std::unordered_map<int, std::vector<GhostElementMetadata>> m_ghostlayerMetadata;

  /** Has a plus fault side */
  bool m_hasPlusFault{false};

  MeshReader(int rank);

  public:
  virtual ~MeshReader();

  const std::vector<Element>& getElements() const;
  const std::vector<Vertex>& getVertices() const;
  const std::map<int, MPINeighbor>& getMPINeighbors() const;
  const std::map<int, std::vector<MPINeighborElement>>& getMPIFaultNeighbors() const;
  const std::unordered_map<int, std::vector<GhostElementMetadata>>& getGhostlayerMetadata() const;
  const std::vector<Fault>& getFault() const;
  bool hasFault() const;
  bool hasPlusFault() const;

  virtual bool inlineTimestepCompute() const { return false; }
  virtual bool inlineClusterCompute() const { return false; }

  void displaceMesh(const Eigen::Vector3d& displacement);

  void computeTimestepIfNecessary(const seissol::SeisSol& seissolInstance);

  // scalingMatrix is stored column-major, i.e.
  // scalingMatrix_ij = scalingMatrix[j][i]
  void scaleMesh(const Eigen::Matrix3d& scalingMatrix);

  /**
   * Reconstruct the fault information from the boundary conditions
   */
  void extractFaultInformation(const VrtxCoords& refPoint,
                               seissol::initializer::parameters::RefPointMethod refPointMethod);

  void exchangeGhostlayerMetadata();
};

} // namespace seissol::geometry

#endif // SEISSOL_SRC_GEOMETRY_MESHREADER_H_

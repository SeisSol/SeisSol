/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 *
 * @section DESCRIPTION
 * Read a mesh file in memory efficient way
 **/

#ifndef MESH_READER_H
#define MESH_READER_H

#include "MeshDefinition.h"

#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "Initializer/Parameters/DRParameters.h"

namespace seissol::geometry {

struct GhostElementMetadata {
  double vertices[4][3];
  int group;
  GlobalElemId globalId;
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

  void displaceMesh(const Eigen::Vector3d& displacement);

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

#endif // MESH_READER_H

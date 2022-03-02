/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
 * Read Gambit Mesh and Metis Partition in memory efficient way
 **/

#ifndef MESH_READER_H
#define MESH_READER_H

#include "MeshDefinition.h"
#include "MeshTools.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

class MeshReader
{
protected:
	const int m_rank;

	std::vector<Element> m_elements;

	std::vector<Vertex> m_vertices;

	/** Convert global element index to local */
	std::map<int, int> m_g2lElements;

	/** Convert global vertex index to local */
	std::map<int, int> m_g2lVertices;

	/** Number of MPI neighbors */
	std::map<int, MPINeighbor> m_MPINeighbors;

	/** Number of MPI fault neighbors */
	std::map<int, std::vector<MPINeighborElement> > m_MPIFaultNeighbors;

	/** Fault information */
	std::vector<Fault> m_fault;

	/** Has a plus fault side */
	bool m_hasPlusFault;

protected:
	MeshReader(int rank)
		: m_rank(rank), m_hasPlusFault(false)
	{}

public:
	virtual ~MeshReader()
	{
	}

	const std::vector<Element>& getElements() const
	{
		return m_elements;
	}

	const std::vector<Vertex>& getVertices() const
	{
		return m_vertices;
	}

	const std::map<int, MPINeighbor>& getMPINeighbors() const
	{
		return m_MPINeighbors;
	}

	const std::map<int, std::vector<MPINeighborElement> >& getMPIFaultNeighbors() const
	{
		return m_MPIFaultNeighbors;
	}

	const std::vector<Fault>& getFault() const
	{
		return m_fault;
	}

	bool hasFault() const
	{
		return m_fault.size() > 0;
	}

	bool hasPlusFault() const
	{
		return m_hasPlusFault;
	}

  void displaceMesh(double const displacement[3])
  {
    for (unsigned vertexNo = 0; vertexNo < m_vertices.size(); ++vertexNo) {
      for (unsigned i = 0; i < 3; ++i) {
        m_vertices[vertexNo].coords[i] += displacement[i];
      }
    }
  }

  // scalingMatrix is stored column-major, i.e.
  // scalingMatrix_ij = scalingMatrix[j][i]
  void scaleMesh(double const scalingMatrix[3][3])
  {
    for (unsigned vertexNo = 0; vertexNo < m_vertices.size(); ++vertexNo) {
      double x = m_vertices[vertexNo].coords[0];
      double y = m_vertices[vertexNo].coords[1];
      double z = m_vertices[vertexNo].coords[2];
      for (unsigned i = 0; i < 3; ++i) {
        m_vertices[vertexNo].coords[i] = scalingMatrix[0][i] * x + scalingMatrix[1][i] * y + scalingMatrix[2][i] * z;
      }
    }
  }

	/**
	 * Reconstruct the fault information from the boundary conditions
	 */
	void findFault(const VrtxCoords refPoint, const int refPointMethod)
	{
		for (std::vector<Element>::iterator i = m_elements.begin();
				i != m_elements.end(); i++) {

			for (int j = 0; j < 4; j++) {
				// Set default mpi fault indices
				i->mpiFaultIndices[j] = -1;

				if (i->boundaries[j] != 3)
					continue;

				// DR boundary

				if (i->neighborRanks[j] == m_rank) {
					// Completely local DR boundary

					if (i->neighbors[j] < i->localId)
						// This was already handled by the other side
						continue;
				} else {
					// Handle boundary faces

					// FIXME we use the MPI number here for the neighbor element id
					// It is not very nice but should generate the correct ordering.
					MPINeighborElement neighbor = {i->localId, j, i->mpiIndices[j], i->neighborSides[j]};
					m_MPIFaultNeighbors[i->neighborRanks[j]].push_back(neighbor);
				}

				Fault f;

				// Detect +/- side
				// Computes the distance between the bary center of the tetrahedron and the face
				// Does not work for all meshes
//				VrtxCoords elementCenter;
//				VrtxCoords faceCenter;
//				MeshTools::center(*i, m_vertices, elementCenter);
//				MeshTools::center(*i, j, m_vertices, faceCenter);
//
//				bool isPlus = (MeshTools::distance(elementCenter, center)
//						< MeshTools::distance(faceCenter, center));

				// Compute normal of the DR face
				// Boundary side vector pointing in chi- and tau-direction
				VrtxCoords chiVec, tauVec;
				MeshTools::sub(m_vertices[i->vertices[MeshTools::FACE2NODES[j][1]]].coords,
						m_vertices[i->vertices[MeshTools::FACE2NODES[j][0]]].coords,
						chiVec);
				MeshTools::sub(m_vertices[i->vertices[MeshTools::FACE2NODES[j][2]]].coords,
						m_vertices[i->vertices[MeshTools::FACE2NODES[j][0]]].coords,
						tauVec);
				MeshTools::cross(chiVec, tauVec, f.normal);

				// Normalize normal
				MeshTools::mul(f.normal, 1.0 / MeshTools::norm(f.normal), f.normal);

				// Check whether the tetrahedron and the reference point are on the same side of the face
				VrtxCoords tmp1, tmp2;
				MeshTools::sub(refPoint, m_vertices[i->vertices[MeshTools::FACE2NODES[j][0]]].coords, tmp1);
				MeshTools::sub(m_vertices[i->vertices[MeshTools::FACE2MISSINGNODE[j]]].coords,
						m_vertices[i->vertices[MeshTools::FACE2NODES[j][0]]].coords,
						tmp2);
				bool isPlus;
				if (refPointMethod == 0) {
					isPlus = MeshTools::dot(tmp1, f.normal) * MeshTools::dot(tmp2, f.normal) > 0;
				} else {
					isPlus = MeshTools::dot(refPoint, f.normal)> 0;
				}
				// Fix normal direction and get correct chiVec
				if (!isPlus) {
					// In case of a minus side, compute chi using node 0 and 1 from the plus side
					MeshTools::sub(m_vertices[i->vertices[MeshTools::FACE2NODES[j]
					              [MeshTools::NEIGHBORFACENODE2LOCAL[(3+1-i->sideOrientations[j])%3]]]].coords,
							m_vertices[i->vertices[MeshTools::FACE2NODES[j]
							      [MeshTools::NEIGHBORFACENODE2LOCAL[(3+0-i->sideOrientations[j])%3]]]].coords,
							chiVec);

          MeshTools::sub(m_vertices[i->vertices[MeshTools::FACE2NODES[j]
					              [MeshTools::NEIGHBORFACENODE2LOCAL[(3+2-i->sideOrientations[j])%3]]]].coords,
							m_vertices[i->vertices[MeshTools::FACE2NODES[j]
							      [MeshTools::NEIGHBORFACENODE2LOCAL[(3+0-i->sideOrientations[j])%3]]]].coords,
							tauVec);

          MeshTools::cross(chiVec, tauVec, f.normal);
          MeshTools::mul(f.normal, 1.0 / MeshTools::norm(f.normal), f.normal);
				}

				// Compute vector inside the triangle's plane for the rotation matrix
				MeshTools::mul(chiVec, 1.0 / MeshTools::norm(chiVec), f.tangent1);
				// Compute second vector in the plane, orthogonal to the normal and tangent 1 vectors
				MeshTools::cross(f.normal, f.tangent1, f.tangent2);

				// Index of the element on the other side
				int neighborIndex = i->neighbors[j] == static_cast<int>(m_elements.size()) ? -1 : i->neighbors[j];

				if (isPlus) {
					f.element = i->localId;
					f.side = j;
					f.neighborElement = neighborIndex;
					f.neighborSide = i->neighborSides[j];
				} else {
					f.element = neighborIndex;
					f.side = i->neighborSides[j];
					f.neighborElement = i->localId;
					f.neighborSide = j;
				}

				m_fault.push_back(f);

				// Check if we have a plus fault side
				if (isPlus || neighborIndex >= 0)
					m_hasPlusFault = true;
			}
		}

		// Sort fault neighbor lists and update MPI fault indices
		for (std::map<int, std::vector<MPINeighborElement> >::iterator i = m_MPIFaultNeighbors.begin();
				i != m_MPIFaultNeighbors.end(); i++) {

			if (i->first > m_rank)
				std::sort(i->second.begin(), i->second.end(), compareLocalMPINeighbor);
			else
				std::sort(i->second.begin(), i->second.end(), compareRemoteMPINeighbor);

			// Set the MPI fault number of all elements
			for (int j = 0; j < static_cast<int>(i->second.size()); j++) {
				m_elements[i->second[j].localElement].mpiFaultIndices[i->second[j].localSide] = j;
			}
		}
	}

protected:
	static bool compareLocalMPINeighbor(const MPINeighborElement &elem1, const MPINeighborElement &elem2)
	{
		return (elem1.localElement < elem2.localElement)
				|| (elem1.localElement == elem2.localElement && elem1.localSide < elem2.localSide);
	}

	static bool compareRemoteMPINeighbor(const MPINeighborElement &elem1, const MPINeighborElement &elem2)
	{
		return (elem1.neighborElement < elem2.neighborElement)
					|| (elem1.neighborElement == elem2.neighborElement && elem1.neighborSide < elem2.neighborSide);
	}
};

#endif // MESH_READER_H

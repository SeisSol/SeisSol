/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#ifndef PUMLREADER_H
#define PUMLREADER_H

#include "MeshReader.h"
#include "Parallel/MPI.h"


// TODO(Lukas) Remove or refactor
#include "PUML/PUML.h"

// from https://stackoverflow.com/a/33965522
#include <array>
void hash_combine(std::size_t& seed, std::size_t value);
using arr_t = std::array<double, 3>;

struct container_hasher {
     template<class T>
     std::size_t operator()(const T& c) const {
         std::size_t seed = 0;
         for(const auto& elem : c) {
             hash_combine(seed, std::hash<typename T::value_type>()(elem));
         }
         return seed;
     }
};

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
      std::size_t h1 = std::hash<T1>{}(p.first);
      std::size_t h2 = std::hash<T2>{}(p.second);

      // TODO(Lukas): Use hash_combine.
      //std::size_t seed = 0;
      //hash_combine(seed, h1);
      //hash_combine(seed, h2);
      //return seed;
      return h1 ^ h2;
    }
};

#ifndef PUML_PUML_H
namespace PUML
{
	class TETPUML;
}
#endif // PUML_PUML_H

namespace seissol {
  namespace initializers {
    namespace time_stepping {
      class LtsWeights;
    }
  }
}

namespace seissol
{
class PUMLReader : public MeshReader
{
public:
        PUMLReader(const char* meshFile, initializers::time_stepping::LtsWeights* ltsWeights = nullptr, double tpwgt = 1.0, bool readPartitionFromFile = false);

private:
	/**
	 * Read the mesh
	 */
	void read(PUML::TETPUML &puml, const char* meshFile);

	/**
	 * Create the partitioning
	 */
	void partition(PUML::TETPUML &puml, initializers::time_stepping::LtsWeights* ltsWeights, double tpwgt, const char *meshFile, bool readPartitionFromFile);
	int readPartition(PUML::TETPUML &puml, int* partition, const char *meshFile);
	void writePartition(PUML::TETPUML &puml, int* partition, const char *meshFile);
	/**
	 * Generate the PUML data structure
	 */
	void generatePUML(PUML::TETPUML &puml);

	/**
	 * Get the mesh
	 */
	void getMesh(const PUML::TETPUML &puml);

	void addMPINeighor(const PUML::TETPUML &puml, int rank, const std::vector<unsigned int> &faces);

private:
	static int FACE_PUML2SEISSOL[4];
	static int FACEVERTEX2ORIENTATION[4][4];
	static int FIRST_FACE_VERTEX[4];


	void getFaceNodes(const Element &element,
			  const std::vector<PUML::TETPUML::vertex_t> &vertices,
			  size_t face,
			  std::array<const double*, 3> &faceNodesCoords
			  ) {
	  auto &faceNodesIds = MeshTools::FACE2NODES[FACE_PUML2SEISSOL[face]];
	  for (int k = 0; k < faceNodesCoords.size(); ++k) {
	    const auto curVertexId = element.vertices[faceNodesIds[k]];
	    faceNodesCoords[k] = vertices[curVertexId].coordinate();
	  }
	}


	  Plane computePlane(const std::array<const double*, 3> &faceNodesCoords) const {
	    Plane plane;
	    // Compute an orthonormal coordinate system for each plane.
	    // Basis is [normal, tangent1, tangent2].
	    VrtxCoords ac;
	    MeshTools::sub(faceNodesCoords[1], faceNodesCoords[0], plane.tangent1);
	    MeshTools::sub(faceNodesCoords[2], faceNodesCoords[0], ac);
	    MeshTools::cross(plane.tangent1, ac, plane.normal);
	    // normalize normal
	    MeshTools::mul(plane.normal,
			   1/(MeshTools::norm(plane.normal)),
			   plane.normal);
	    // take absolut value of normal to standarize direction
	    // This way, opposite planes have the same coordinate system.
	    for (int k = 0; k < 3; ++k) {
	      plane.normal[k] = std::abs(plane.normal[k]);
	    }
	    if (std::abs(plane.normal[0]) > std::abs(plane.normal[1])) {
	      plane.tangent1[0] = plane.normal[2];
	      plane.tangent1[1] = 0.0;
	      plane.tangent1[2] = -plane.normal[0];
	    } else {
	      plane.tangent1[0] = 0.0;
	      plane.tangent1[1] = -plane.normal[2];
	      plane.tangent1[2] = plane.normal[1];
	    }
	    MeshTools::cross(plane.normal, plane.tangent1, plane.tangent2);

	    for (int k = 0; k < 3; ++k) {
	      // A point on the plane is needed for the plane equations.
	      // This also differentiates opposite planes.
	      plane.point[k] = faceNodesCoords[0][k];
	    }
	    return plane;
	  }

	  int findMatchingCell(int cellId,
			       int faceId,
			       Plane facePlane,
			       std::array<const double*, 3> &faceNodesCoords,
			       const std::vector<Plane> &planes,
			       const std::vector<std::unordered_set<PeriodicVertex>> &planes2Vertices,
			       std::unordered_map<size_t, std::unordered_set<int>>&vertices2Cells,
			       std::unordered_map<size_t, std::unordered_set<size_t>> &vertices2Vertices) {
	    auto &faceNodesIds = MeshTools::FACE2NODES[FACE_PUML2SEISSOL[faceId]];
	    auto result = -1;

	    auto maxCounterAllPlanes = -1; // TODO(Lukas) Name
	    
	    bool foundPlane = false;
	    for (int k = 0; k < planes.size(); ++k) {
	      const auto &curPlane = planes[k];
	      if (curPlane == facePlane) continue;

	      // Move current plane s.t. it cuts the other plane at its point.
	      // Remember: opposite planes have different offsets but same coordinate system!
	      //if (cellId < 0) logInfo(0) << "Checking plane" << k ;
	      auto movedPlane = planes[k];
	      std::copy_n(facePlane.point, 3, movedPlane.point);
	      //if (cellId < 0) logInfo(0) << "init moved plane" ;
	      bool containsAllPoints = true;
	      if (cellId < 0) logInfo(0) << "Checking planes";
	      for (int l = 0; l < 3; ++l) {
		containsAllPoints &=
		  movedPlane.containsPoint(faceNodesCoords[l]);
	      }
	      //if (cellId < 0) logInfo(0) << "Contains all points =" << containsAllPoints;
	      if (containsAllPoints) {
		foundPlane = true;
		// We count to which face most points belong.
		auto faceCounter = std::unordered_map<size_t, int>{};
		// Iterate over all points of the edge.
		for (int l = 0; l < 3; ++l) {
		  const auto* curCoords = faceNodesCoords[l];
		  const auto curVertex = PeriodicVertex{faceNodesCoords[l],
							planes[k]};
		  const auto curVertexId = (cellId > 0) ? m_elements[cellId].vertices[faceNodesIds[l]] : -1;
		  // We compare with all other points in the plane.
		  // Note that this is O(N^2); this doesn't really matter,
		  // because normally N is relatively small for a single node.
		  for (const auto &other : planes2Vertices[k]) {
		    /*
		    std::cout << "Checking vertex: "
		      << "x = " << curCoords[0] << "|" << other.localCoord[0]
		      << "\t y =" << curCoords[1] << "|" << other.localCoord[1]
                      << "\t z =" << curCoords[2] << "|" << other.localCoord[2]
			<< std::endl;
		    */

		    if (curVertex.isApproxEq(other)) {
		      //std::cout << "Found matching vertex." << std::endl;
		      vertices2Vertices[curVertexId].insert(other.vertexId); // TODO not if -1
		      assert(curVertexId < 0 || vertices2Cells.find(curVertexId) !=
			     vertices2Cells.end());
		      assert(vertices2Cells.find(other.vertexId) != vertices2Cells.end());
		      for (const auto otherCellId : vertices2Cells[other.vertexId]) {
			// std::cout << "Found matching cell." << std::endl;
			faceCounter[otherCellId] += 1;
		      }
		    }
		  }
		}

		// Check which cell was mentioned most often
		// This is the matching cell, as most vertices belong to it.
		auto bestCellId = -1;
		auto maxCounter = -1;
		for (const auto& it : faceCounter) {
		  if (it.first != cellId && it.second > maxCounter) {
		    maxCounter = it.second;
		    bestCellId = it.first;
		  }
		}
		std::cout << "Result = " << result
			  << " bestCellId = " << bestCellId
			  << " maxCounter = " << maxCounter
			  << " maxCounterAllPlanes = " << maxCounterAllPlanes
			  << std::endl;
		// Make sure that we have an unique matching.
		assert(!(maxCounter >= 3 && maxCounterAllPlanes >= 3));
		if (maxCounter > maxCounterAllPlanes) {
		  result = bestCellId;
		}
		maxCounterAllPlanes = std::max(maxCounter, maxCounterAllPlanes);

		// For MPI, it is often needed to try out multiple planes.
		// TODO(Lukas) Why?
		//result = std::max(result, bestCellId);
		if (faceCounter.size() == 0) {
		  std::cout << "No face mentioned." << std::endl;
		}
		std::cout << "Result = " << result << std::endl;
		/*
		std::cout << "Return matching " << result
			  << " FaceCounter size " << faceCounter.size()
			  << std::endl;
		*/
		//return result; // TODO(Lukas) Is this correct?
	      }
	    }
	    if (!foundPlane) {
	      std::cout << "Did not find plane!" << std::endl;
	    }
	    std::cout << "\n" << std::endl;
	    return result;
	  }
};
 
}

#endif // PUMLREADER_H

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
#include "PUML/PUML.h"

namespace seissol {
  namespace initializers {
    namespace time_stepping {
      class LtsWeights;
    }
  }
}

namespace seissol::geometry
{
class PUMLReader : public seissol::geometry::MeshReader
{
public:
        PUMLReader(const char* meshFile, double maximumAllowedTimeStep, const char* checkPointFile,
            initializers::time_stepping::LtsWeights* ltsWeights = nullptr, double tpwgt = 1.0, bool readPartitionFromFile = false);

private:
	/**
	 * Read the mesh
	 */
	void read(PUML::TETPUML &puml, const char* meshFile);

	/**
	 * Create the partitioning
	 */
	void partition(PUML::TETPUML &puml, initializers::time_stepping::LtsWeights* ltsWeights, double tpwgt, const char *meshFile, bool readPartitionFromFile, const char* checkPointFile);
	int readPartition(PUML::TETPUML &puml, int* partition, const char *checkPointFile);
	void writePartition(PUML::TETPUML &puml, int* partition, const char *checkPointFile);
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
};

/*
 * Possible types of boundary conditions for SeisSol.
 */
enum class BCType {
	internal, external, unknown
};

/**
 * Decodes the boundary condition tag into a BCType.
 */
constexpr BCType bcToType(int id) {
	if (id == 0 || id == 3 || id > 64) {
		return BCType::internal;
	} else if (id == 1 || id == 2 || id == 4 || id == 5 || id == 6 || id == 7) {
		return BCType::external;
	} else {
		return BCType::unknown;
	}
}

/**
 * Decodes the boundary condition tag into a string representation.
 */
inline std::string bcToString(int id) {
	if (id == 0) { return std::string("regular"); }
	else if (id == 1) { return std::string("free surface"); }
	else if (id == 2) { return std::string("free surface with gravity"); }
	else if (id == 3) { return std::string("dynamic rupture"); }
        else if (id == 4) { return std::string("dirichlet"); }
	else if (id == 5) { return std::string("absorbing"); }
	else if (id == 6) { return std::string("periodic"); }
        else if (id == 7) { return std::string("analytic"); }
	else if (id > 64) {
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
inline bool checkMeshCorrectnessLocally(PUML::TETPUML::face_t face, int* cellNeighbors, int side, int sideBC, int cellIdAsInFile) {
	// if a face is an internal face, it has to have a neighbor on either this rank or somewhere else:
	if (bcToType(sideBC) == BCType::internal) {
		if (cellNeighbors[side] < 0 && !face.isShared()) {
			logWarning() << "Element" << cellIdAsInFile << ", side" << side << " has a"
				<< bcToString(sideBC)
				<< "boundary condition, but the neighboring element doesn't exist";
			return false;
		}
	}
	// external boundaries must not have neighboring elements:
	else if (bcToType(sideBC) == BCType::external) {
		if (cellNeighbors[side] >= 0 || face.isShared()) {
			logWarning() << "Element" << cellIdAsInFile << ", side" << side << " has a"
				<< bcToString(sideBC)
				<< "boundary condition, but the neighboring element is not flagged -1";
		return false;
	}
	}
	// ignore unknown boundary conditions and warn
	else {
		logWarning() << "Element" << cellIdAsInFile << ", side" << side
			<< " has a boundary condition, which I don't understand" << sideBC;
	  return true;
	}
	return true;
}

}

#endif // PUMLREADER_H

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2017, SeisSol Group
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


#include <algorithm>
#include <cstring>
#include <iostream>

#include "utils/logger.h"

#include "SeisSol.h"
#include "MeshReaderFBinding.h"
#include "GambitReader.h"
#ifdef USE_NETCDF
#include "NetcdfReader.h"
#endif // USE_NETCDF
#if defined(USE_METIS) && defined(USE_HDF) && defined(USE_MPI)
#include "PUMLReader.h"
#endif // defined(USE_METIS) && defined(USE_HDF) && defined(USE_MPI)
#include "Modules/Modules.h"
#include "Monitoring/instrumentation.fpp"
#include "Monitoring/Stopwatch.h"
#include "Numerical_aux/Statistics.h"
#include "Initializer/time_stepping/LtsWeights.h"
#include "Solver/time_stepping/MiniSeisSol.h"

void read_mesh(int rank, MeshReader &meshReader, bool hasFault, double const displacement[3], double const scalingMatrix[3][3])
{
	logInfo(rank) << "Reading mesh. Done.";

	meshReader.displaceMesh(displacement);
	meshReader.scaleMesh(scalingMatrix);

	const std::vector<Element>& elements = meshReader.getElements();
	const std::vector<Vertex>& vertices = meshReader.getVertices();
	const std::map<int, MPINeighbor>& mpiNeighbors = meshReader.getMPINeighbors();

	// Compute maximum element for one vertex
	size_t maxElements = 0;
	for (std::vector<Vertex>::const_iterator i = vertices.begin();
			i != vertices.end(); i++)
		maxElements = std::max(maxElements, i->elements.size());

	allocelements(elements.size());
	allocvertices(vertices.size(), maxElements);
	allocbndobjs(mpiNeighbors.size());

	// Set vertices
	int size;
	double* verticesXY;
	int* verticesNElements;
	int* verticesElements;
	getverticesxy(&size, &verticesXY);
	getverticesnelements(&size, &verticesNElements);
	getverticeselements(&size, &verticesElements);
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < 3; j++) {
			verticesXY[i*3+j] = vertices[i].coords[j];
		}

		verticesNElements[i] = vertices[i].elements.size();

		for (unsigned int j = 0; j < vertices[i].elements.size(); j++) {
			verticesElements[i+j*vertices.size()] = vertices[i].elements[j] + 1;
		}
	}

	// Set elements
	int* elementVertices;
	int* sideNeighbor;
	int* localNeighborSide;
	int* localNeighborVrtx;
	int* reference;
	int* mpiReference;
	int* mpiNumber;
	int* boundaryToObject;
	getelementvertices(&size, &elementVertices);
	getsideneighbor(&size, &sideNeighbor);
	getlocalneighborside(&size, &localNeighborSide);
	getlocalneighborvrtx(&size, &localNeighborVrtx);
	getreference(&size, &reference);
	getmpireference(&size, &mpiReference);
	getmpinumber(&size, &mpiNumber);
	getboundarytoobject(&size, &boundaryToObject);

	for (int i = 0; i < size; i++) {
		reference[i*5] = elements[i].material;

		for (int j = 0; j < 4; j++) {
			elementVertices[i*4+j] = elements[i].vertices[j] + 1;

			sideNeighbor[i*4+j] = elements[i].neighbors[j] + 1;

			reference[i*5+j+1] = elements[i].boundaries[j];
			switch(reference[i*5+j+1]) {
			case 0:
			case 3:
			case 6:
				localNeighborSide[i*4+j] = elements[i].neighborSides[j] + 1;
				localNeighborVrtx[i*4+j] = elements[i].sideOrientations[j] + 1;
				boundaryToObject[i*4+j] = 0;
				break;
			default:
				localNeighborSide[i*4+j] = 1;
				localNeighborVrtx[i*4+j] = 1;
				boundaryToObject[i*4+j] = 1;
			}

			mpiReference[i*5+j+1] = (elements[i].neighborRanks[j] == rank ? 0 : 1);
			if (mpiReference[i*5+j+1]) {
				boundaryToObject[i*4+j] = mpiNeighbors.at(elements[i].neighborRanks[j]).localID + 1;
				mpiNumber[i*4+j] = elements[i].mpiIndices[j] + 1;
			} else {
				mpiNumber[i*4+j] = -1;
			}
		}
	}

	// Set Boundary objects
	for (std::map<int, MPINeighbor>::const_iterator i = mpiNeighbors.begin(); i != mpiNeighbors.end(); i++) {
		int localID = i->second.localID;
		int* bnddomainelements;
		setbndnelem(localID + 1, i->second.elements.size());

		allocbndobj(localID + 1, i->second.elements.size());

		setbndrank(localID + 1, i->first);
		getbnddomainelements(localID + 1, &size, &bnddomainelements);

		for (int j = 0; j < size; j++) {
			bnddomainelements[j] = i->second.elements[j].localElement + 1;
		}
	}

	if (hasFault) {
		logInfo(rank) << "Extracting fault information";

		VrtxCoords center;
		int refPointMethod;
		getfaultreferencepoint(&center[0], &center[1], &center[2], &refPointMethod);
		meshReader.findFault(center, refPointMethod);

		int* mpiNumberDr;
		getmpinumberdr(&size, &mpiNumberDr);

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < 4; j++) {
				if (elements[i].neighborRanks[j] == rank)
					mpiNumberDr[i*4+j] = -1;
				else
					mpiNumberDr[i*4+j] = elements[i].mpiFaultIndices[j] + 1;
			}
		}

		const std::vector<Fault>& fault = meshReader.getFault();
		const std::map<int, std::vector<MPINeighborElement> >& mpiFaultNeighbors = meshReader.getMPIFaultNeighbors();

		allocfault(fault.size());

		if (meshReader.hasPlusFault())
			hasplusfault();

		if (fault.size() > 0) {
			int* faultface;
			double* faultnormals;
			double* faulttangent1;
			double* faulttangent2;
			getfaultface(&size, &faultface);
			getfaultnormals(&size, &faultnormals);
			getfaulttangent1(&size, &faulttangent1);
			getfaulttangent2(&size, &faulttangent2);

			for (int i = 0; i < size; i++) {
				faultface[i] = fault[i].element + 1;
				faultface[i + size] = fault[i].side + 1;
				faultface[i + size*2] = fault[i].neighborElement + 1;
				faultface[i + size*3] = fault[i].neighborSide + 1;

				memcpy(&faultnormals[i*3], fault[i].normal, sizeof(double)*3);
				memcpy(&faulttangent1[i*3], fault[i].tangent1, sizeof(double)*3);
				memcpy(&faulttangent2[i*3], fault[i].tangent2, sizeof(double)*3);
			}
		}

		for (std::map<int, std::vector<MPINeighborElement> >::const_iterator i = mpiFaultNeighbors.begin();
				i != mpiFaultNeighbors.end(); i++) {
			int localID = mpiNeighbors.at(i->first).localID;

			int* bndfaultelements;
			setbndfaultnelem(localID + 1, i->second.size());

			allocbndobjfault(localID + 1, i->second.size());

			getbndfaultelements(localID + 1, &size, &bndfaultelements);

			for (int j = 0; j < size; j++) {
				bndfaultelements[j] = i->second[j].localElement + 1;
			}
		}
	}

	seissol::SeisSol::main.getLtsLayout().setMesh(meshReader);

	// Setup the communicator for dynamic rupture
	seissol::MPI::mpi.fault.init(meshReader.getFault().size() > 0);

	// Call the post mesh initialization hook
	seissol::Modules::callHook<seissol::POST_MESH>();
}

extern "C" {

void read_mesh_gambitfast_c(int rank, const char* meshfile, const char* partitionfile, bool hasFault, double const displacement[3], double const scalingMatrix[3][3])
{
	SCOREP_USER_REGION("read_mesh", SCOREP_USER_REGION_TYPE_FUNCTION);

	logInfo(rank) << "Reading Gambit mesh using fast reader";
	logInfo(rank) << "Parsing mesh and partition file:" << meshfile << ';' << partitionfile;

	Stopwatch watch;
	watch.start();

	seissol::SeisSol::main.setMeshReader(new GambitReader(rank, meshfile, partitionfile));

	read_mesh(rank, seissol::SeisSol::main.meshReader(), hasFault, displacement, scalingMatrix);

	watch.pause();
	watch.printTime("Mesh initialized in:");
}

void read_mesh_netcdf_c(int rank, int nProcs, const char* meshfile, bool hasFault, double const displacement[3], double const scalingMatrix[3][3])
{
	SCOREP_USER_REGION("read_mesh", SCOREP_USER_REGION_TYPE_FUNCTION);

#ifdef USE_NETCDF
	logInfo(rank) << "Reading netCDF mesh" << meshfile;

	Stopwatch watch;
	watch.start();

	seissol::SeisSol::main.setMeshReader(new NetcdfReader(rank, nProcs, meshfile));

	read_mesh(rank, seissol::SeisSol::main.meshReader(), hasFault, displacement, scalingMatrix);

	watch.pause();
	watch.printTime("Mesh initialized in:");

#else // USE_NETCDF
	logError() << "netCDF not supported";
#endif // USE_NETCDF
}


void read_mesh_puml_c(const char* meshfile, const char* checkPointFile, bool hasFault, double const displacement[3], double const scalingMatrix[3][3], char const* easiVelocityModel, int clusterRate)
{
	SCOREP_USER_REGION("read_mesh", SCOREP_USER_REGION_TYPE_FUNCTION);

#if defined(USE_METIS) && defined(USE_HDF) && defined(USE_MPI)
	const int rank = seissol::MPI::mpi.rank();
  	double tpwgt = 1.0;
	if (seissol::MPI::mpi.size() > 1) {
	  logInfo(rank) << "Running mini SeisSol to determine node weight";
	  tpwgt = 1.0 / seissol::miniSeisSol(seissol::SeisSol::main.getMemoryManager());

	  const auto summary = seissol::statistics::parallelSummary(tpwgt);
	  logInfo(rank) << "Node weights: mean =" << summary.mean
			<< " std =" << summary.std
			<< " min =" << summary.min
			<< " median =" << summary.median
			<< " max =" << summary.max;
	}
	
	logInfo(rank) << "Reading PUML mesh" << meshfile;

	Stopwatch watch;
	watch.start();

	bool readPartitionFromFile = seissol::SeisSol::main.simulator().checkPointingEnabled();

	seissol::initializers::time_stepping::LtsWeights ltsWeights(easiVelocityModel, clusterRate);
	seissol::SeisSol::main.setMeshReader(new seissol::PUMLReader(meshfile, checkPointFile, &ltsWeights, tpwgt, readPartitionFromFile));

	read_mesh(rank, seissol::SeisSol::main.meshReader(), hasFault, displacement, scalingMatrix);

	watch.pause();
	watch.printTime("Mesh initialized in:");

#else // defined(USE_METIS) && defined(USE_HDF) && defined(USE_MPI)
	logError() << "PUML is currently only supported for MPI";
#endif // defined(USE_METIS) && defined(USE_HDF) && defined(USE_MPI)
}

}

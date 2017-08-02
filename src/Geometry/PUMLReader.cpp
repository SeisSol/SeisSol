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

#include <algorithm>
#include <cassert>
#include <string>
#include <unordered_map>

#include "PUML/PUML.h"
#include "PUML/PartitionMetis.h"
#include "PUML/Downward.h"
#include "PUML/Neighbor.h"

#include "PUMLReader.h"
#include "Monitoring/instrumentation.fpp"

class GlobalFaceSorter
{
private:
	const PUML::TETPUML &m_puml;

public:
	GlobalFaceSorter(const PUML::TETPUML &puml)
		: m_puml(puml)
	{ }

	bool operator()(unsigned int a, unsigned int b) const
	{
		return m_puml.faces()[a].gid() < m_puml.faces()[b].gid();
	}
};

/**
 * @todo Cleanup this code
 */
seissol::PUMLReader::PUMLReader(const char *meshFile)
	: MeshReader(MPI::mpi.rank())
{
	PUML::TETPUML puml;
	puml.setComm(MPI::mpi.comm());

	read(puml, meshFile);

	partition(puml);

	generatePUML(puml);

	getMesh(puml);
}

void seissol::PUMLReader::read(PUML::TETPUML &puml, const char* meshFile)
{
	SCOREP_USER_REGION("PUMLReader_read", SCOREP_USER_REGION_TYPE_FUNCTION);

	std::string file(meshFile);

	puml.open((file + ":/connect").c_str(), (file + ":/geometry").c_str());
	puml.addData((file + ":/group").c_str(), PUML::CELL);
	puml.addData((file + ":/boundary").c_str(), PUML::CELL);
}

void seissol::PUMLReader::partition(PUML::TETPUML &puml)
{
	SCOREP_USER_REGION("PUMLReader_partition", SCOREP_USER_REGION_TYPE_FUNCTION);

	PUML::TETPartitionMetis metis(puml.originalCells(), puml.numOriginalCells());
	int* partition = new int[puml.numOriginalCells()];
	metis.partition(partition);

	puml.partition(partition);
	delete [] partition;
}

void seissol::PUMLReader::generatePUML(PUML::TETPUML &puml)
{
	SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

	puml.generateMesh();
}

void seissol::PUMLReader::getMesh(const PUML::TETPUML &puml)
{
	SCOREP_USER_REGION("PUMLReader_getmesh", SCOREP_USER_REGION_TYPE_FUNCTION);

	const int rank = MPI::mpi.rank();

	const std::vector<PUML::TETPUML::cell_t> &cells = puml.cells();
	const std::vector<PUML::TETPUML::face_t> &faces = puml.faces();
	const std::vector<PUML::TETPUML::vertex_t> &vertices = puml.vertices();

	const int* material = puml.cellData(0);
	const int* boundaryCond = puml.cellData(1);

	std::unordered_map<int, std::vector<unsigned int> > neighborInfo; // List of shared local face ids

	// Compute everything local
	m_elements.resize(cells.size());
	for (unsigned int i = 0; i < cells.size(); i++) {
		m_elements[i].localId = i;

		// Vertices
		PUML::Downward::vertices(puml, cells[i], reinterpret_cast<unsigned int*>(m_elements[i].vertices));

		// Neighbor information
		unsigned int faceids[4];
		PUML::Downward::faces(puml, cells[i], faceids);
		int neighbors[4];
		PUML::Neighbor::face(puml, i, neighbors);
		for (unsigned int j = 0; j < 4; j++) {
			if (neighbors[j] < 0) {
				m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = cells.size();

				if (!faces[faceids[j]].isShared()) {
					// Boundary sides
					m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = rank;
				} else {
					// MPI Boundary
					neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

					m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = faces[faceids[j]].shared()[0];
				}
			} else {
				assert(neighbors[j] < cells.size());

				m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = neighbors[j];

				int nfaces[4];
				PUML::Neighbor::face(puml, neighbors[j], nfaces);
				int* back = std::find(nfaces, nfaces+4, i);
				assert(back < nfaces+4);

				m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]] = FACE_PUML2SEISSOL[back-nfaces];

				const unsigned int firstVertex = m_elements[i].vertices[FIRST_FACE_VERTEX[FACE_PUML2SEISSOL[j]]];

				unsigned int nvertices[4];
				PUML::Downward::vertices(puml, cells[neighbors[j]], nvertices);
				unsigned int* neighborFirstVertex = std::find(nvertices, nvertices+4, firstVertex);

				m_elements[i].sideOrientations[FACE_PUML2SEISSOL[j]]
					= FACEVERTEX2ORIENTATION[m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]]][neighborFirstVertex-nvertices];
				assert(m_elements[i].sideOrientations[FACE_PUML2SEISSOL[j]] >= 0);

				m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = rank;
			}

			m_elements[i].boundaries[FACE_PUML2SEISSOL[j]] = (boundaryCond[i] >> (j*8)) & 0xFF;

			m_elements[i].mpiIndices[FACE_PUML2SEISSOL[j]] = 0;
		}

		m_elements[i].material = material[i];
	}

	// Exchange ghost layer information and generate neighbor list
	char** copySide = new char*[neighborInfo.size()];
	char** ghostSide = new char*[neighborInfo.size()];
	unsigned long** copyFirstVertex = new unsigned long*[neighborInfo.size()];
	unsigned long** ghostFirstVertex = new unsigned long*[neighborInfo.size()];

	MPI_Request* requests = new MPI_Request[neighborInfo.size() * 4];

	std::unordered_set<unsigned int> t;
	unsigned int sum = 0;
	unsigned int k = 0;
	for (std::unordered_map<int, std::vector<unsigned int> >::iterator it = neighborInfo.begin();
			it != neighborInfo.end(); ++it, k++) {
		// Need to sort the neighborInfo vectors once
		GlobalFaceSorter faceSorter(puml);
		std::sort(it->second.begin(), it->second.end(), faceSorter);

		t.insert(it->second.begin(), it->second.end());
		sum += it->second.size();

		// Create MPI neighbor list
		addMPINeighor(puml, it->first, it->second);

		copySide[k] = new char[it->second.size()];
		ghostSide[k] = new char[it->second.size()];
		copyFirstVertex[k] = new unsigned long[it->second.size()];
		ghostFirstVertex[k] = new unsigned long[it->second.size()];

		MPI_Irecv(ghostSide[k], it->second.size(), MPI_CHAR, it->first, 0, MPI::mpi.comm(), &requests[k]);
		MPI_Irecv(ghostFirstVertex[k], it->second.size(), MPI_UNSIGNED_LONG, it->first, 0, MPI::mpi.comm(),
			&requests[neighborInfo.size() + k]);

		// Neighbor side
		for (unsigned int i = 0; i < it->second.size(); i++) {
			// The side of boundary
			int cellIds[2];
			PUML::Upward::cells(puml, faces[it->second[i]], cellIds);
			int side = PUML::Downward::faceSide(puml, cells[cellIds[0]], it->second[i]);
			assert(side >= 0 && side < 4);
			copySide[k][i] = side;

			// First vertex of the face on the boundary
			const unsigned int firstVertex
				= m_elements[cellIds[0]].vertices[FIRST_FACE_VERTEX[FACE_PUML2SEISSOL[side]]];
			copyFirstVertex[k][i] = vertices[firstVertex].gid();

			// Set the MPI index
			assert(m_elements[cellIds[0]].mpiIndices[FACE_PUML2SEISSOL[side]] == 0);
			m_elements[cellIds[0]].mpiIndices[FACE_PUML2SEISSOL[side]] = i;
		}

		MPI_Isend(copySide[k], it->second.size(), MPI_CHAR, it->first, 0, MPI::mpi.comm(),
			&requests[neighborInfo.size() * 2 + k]);
		MPI_Isend(copyFirstVertex[k], it->second.size(), MPI_UNSIGNED_LONG, it->first, 0, MPI::mpi.comm(),
			&requests[neighborInfo.size() * 3 + k]);
	}
	assert(t.size() == sum);

	MPI_Waitall(neighborInfo.size() * 4, requests, MPI_STATUSES_IGNORE);

	k = 0;
	for (std::unordered_map<int, std::vector<unsigned int> >::const_iterator it = neighborInfo.begin();
			it != neighborInfo.end(); ++it, k++) {
		for (unsigned int i = 0; i < it->second.size(); i++) {
			// Set neighbor side
			int cellIds[2];
			PUML::Upward::cells(puml, faces[it->second[i]], cellIds);
			assert(cellIds[1] < 0);

			int side = copySide[k][i];

			m_elements[cellIds[0]].neighborSides[FACE_PUML2SEISSOL[side]] = FACE_PUML2SEISSOL[ghostSide[k][i]];

			// Set side sideOrientation
			unsigned long nvertices[4];
			PUML::Downward::gvertices(puml, cells[cellIds[0]], nvertices);

			unsigned long* localFirstVertex = std::find(nvertices, nvertices+4, ghostFirstVertex[k][i]);
			assert(localFirstVertex != nvertices+4);

			m_elements[cellIds[0]].sideOrientations[FACE_PUML2SEISSOL[side]]
				= FACEVERTEX2ORIENTATION[FACE_PUML2SEISSOL[side]][localFirstVertex-nvertices];
			assert(m_elements[cellIds[0]].sideOrientations[FACE_PUML2SEISSOL[side]] >= 0);
		}

		delete [] copySide[k];
		delete [] ghostSide[k];
		delete [] copyFirstVertex[k];
		delete [] ghostFirstVertex[k];
	}

	delete [] copySide;
	delete [] ghostSide;
	delete [] copyFirstVertex;
	delete [] ghostFirstVertex;
	delete [] requests;

	// Set vertices
	m_vertices.resize(vertices.size());
	for (unsigned int i = 0; i < vertices.size(); i++) {
		memcpy(m_vertices[i].coords, vertices[i].coordinate(), 3*sizeof(double));

		PUML::Upward::cells(puml, vertices[i], m_vertices[i].elements);
	}
}

void seissol::PUMLReader::addMPINeighor(const PUML::TETPUML &puml, int rank, const std::vector<unsigned int> &faces)
{
	unsigned int id = m_MPINeighbors.size();
	MPINeighbor& neighbor = m_MPINeighbors[rank];

	neighbor.localID = id;
	neighbor.elements.resize(faces.size());

	for (unsigned int i = 0; i < faces.size(); i++) {
		int cellIds[2];
		PUML::Upward::cells(puml, puml.faces()[faces[i]], cellIds);

		neighbor.elements[i].localElement = cellIds[0];
	}
}


int seissol::PUMLReader::FACE_PUML2SEISSOL[4] = {0, 1, 3, 2};

int seissol::PUMLReader::FACEVERTEX2ORIENTATION[4][4] = {
	{0, 2, 1, -1}, {0, 1, -1, 2}, {0, -1, 2, 1}, {-1, 0, 1, 2}
};

int seissol::PUMLReader::FIRST_FACE_VERTEX[4] = {0, 0, 0, 1};

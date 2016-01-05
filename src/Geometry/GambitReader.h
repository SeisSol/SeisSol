/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013-2014, SeisSol Group
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

#ifndef GAMBIT_READER_H
#define GAMBIT_READER_H

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "MeshReader.h"
#include "SeisSol.h"
#include "vectors.h"

#include "utils/stringutils.h"
#include "utils/logger.h"

class GambitReader : public MeshReader
{
private:
	std::ifstream m_mesh;
#ifdef PARALLEL
	std::ifstream m_partition;
#endif // PARALLEL

	/** Seek positions */
	size_t m_seekVertices;
	size_t m_seekElements;

	int m_nGlobElements;
	int m_nGlobVertices;

public:
	GambitReader(int rank, const char* meshFile, const char* partitionFile)
		: MeshReader(rank), m_mesh(meshFile)
#ifdef PARALLEL
			, m_partition(partitionFile)
#endif // PARALLEL
	{
		// Files are readable?
		if (!m_mesh)
			logError() << "Could not open mesh file" << meshFile;
#ifdef PARALLEL
		if (!m_partition)
			logError() << "Could not open partition file" << partitionFile;
#endif // PARALLEL

		std::string line;

		// Read header information
		getline(m_mesh, line);	// First line contains version
								// we ignore this line for know
		line.clear();
		getline(m_mesh, line);
		utils::StringUtils::trim(line);
		if (line != GAMBIT_FILE_ID)
			logError() << "Not a Gambit mesh file:" << meshFile;

		std::string name;
		getline(m_mesh, name);		// Internal name
		//trim(name);

		getline(m_mesh, line); 	// Skip line:
								// PROGRAM: Gambit VERSION: x.y.z

		getline(m_mesh, line);		// Date
		//trim(line);

		getline(m_mesh, line); 	// Skip problem size names

		int m_nGroups, boundarys, dimensions;
		m_mesh >> m_nGlobVertices;
		m_mesh >> m_nGlobElements;
		m_mesh >> m_nGroups;
		m_mesh >> boundarys;	// Not used at the moment
		m_mesh >> dimensions;
		getline(m_mesh, line); 	// Skip rest of the line
		if (dimensions != 3)
			logError() << "Gambit file does not contain a 3 dimensional mesh";

		getline(m_mesh, line);
		utils::StringUtils::trim(line);
		if (line != ENDSECTION)
			abort();

		// Count size of the local partition
		int l = 0;
		for (int i = 0; i < m_nGlobElements; i++) {
			int elementRank = nextRank();
			if (elementRank == rank)
				m_g2lElements[i] = l++;
		}

		// Find the seek positions for all sections, read local elements and find local vertices
		m_mesh >> std::ws; // Skip white spaces
		while (!m_mesh.eof()) {
			getline(m_mesh, line);

			if (line.find(NODAL_COORDINATES) == 0) {
				// Start coords section
				m_seekVertices = m_mesh.tellg();
			} else if (line.find(ELEMENT_CELLS) == 0) {
				// Start element section
				m_seekElements = m_mesh.tellg();
				parseLocalElements();
			} else if (line.find(ELEMENT_GROUP) == 0) {
				// Start of one group section
				parseLocalGroups();
			} else if (line.find(BOUNDARY_CONDITIONS) == 0) {
				// Start of one boundary section
				parseLocalBoundaries();
			}

			m_mesh >> std::ws; // Skip white spaces
		}

		// Find neighbor MPI ranks
		parseMPINeighborElements();

		// Find local neighbor elements
		parseLocalNeighborElements();

		// Find global to local mapping of vertices
		parseLocalCoordinates();

		// Translate global vertices to local
		translateG2LVertices();

		// Sort neighbor lists
		sortMPINeighborElements();
	}

	virtual ~GambitReader()
	{
	}

private:
	void parseLocalElements()
	{
#ifdef PARALLEL
		m_partition.clear();
		m_partition.seekg(0);
#endif // PARALLEL
		// Do not seek the mesh file, we should be at the correct position

		int numPartitions = seissol::MPI::mpi.size();

		m_elements.resize(m_g2lElements.size());
		int k = 0;
		for (int i = 0; i < m_nGlobElements; i++) {
			std::string line;
			int elementRank = nextRank();

			if (elementRank == m_rank) {
				int n, t;

				assert(static_cast<size_t>(k) < m_elements.size());

				m_elements[k].localId = k;

				m_mesh >> n; // Element number
				m_mesh >> t; // Type
				if (t != 6)
					abort();
				m_mesh >> n; // Number of points
				for (int j = 0; j < 4; j++) {
					m_mesh >> m_elements[k].vertices[j];
					m_elements[k].vertices[j]--;

					if (m_g2lVertices.find(m_elements[k].vertices[j]) == m_g2lVertices.end()) {
						// First time we see this vertex
						// put it in the g2l map and resize the vertex list
						int nextVertex = m_g2lVertices.size();
						m_g2lVertices[m_elements[k].vertices[j]] = nextVertex;

						assert(m_vertices.size() == m_g2lVertices.size() - 1);

						m_vertices.resize(m_g2lVertices.size());
					}

					// Add this element to the vertex list
					m_vertices[m_g2lVertices[m_elements[k].vertices[j]]].elements.push_back(k);

					m_elements[k].rank = m_rank;
					m_elements[k].neighbors[j] = -1;
					m_elements[k].boundaries[j] = 0;
				}

				k++;
				m_mesh >> std::ws; // Skip rest of the line
			} else if (elementRank >= numPartitions) {
				logError() << "Invalid Partition file. Found element rank" << elementRank;
			} else {
				getline(m_mesh, line);
			}

#ifdef PARALLEL
			m_partition >> std::ws;
#endif // PARALLEL
		}

		std::string line;
		m_mesh >> std::ws;
		getline(m_mesh, line);

		if (line.find(ENDSECTION) != 0)
			abort();
	}

	void parseLocalBoundaries()
	{
		static const int faceG2S[4] = {0, 1, 3, 2};

		// Do not seek mesh file, we should already be at the correct position

		int boundaryCondition;
		m_mesh >> boundaryCondition; // Boundary type
		boundaryCondition -= 100;
		int t;
		m_mesh >> t; // Ignore
		int k;
		m_mesh >> k; // Number of boundary conditions

		std::string line;
		getline(m_mesh, line); // Skip rest of the line

		for (int i = 0; i < k; i++) {
			int n;
			m_mesh >> n; // Element number
			n--;

			std::map<int, int>::const_iterator e = m_g2lElements.find(n);

			if (e != m_g2lElements.end()) {
				int t;
				m_mesh >> t; // Element type
				if (t != 6)
					abort();

				int s;
				m_mesh >> s; // Element side
				s = faceG2S[s-1];

				if (boundaryCondition != 3)
					// We still need to find the neighboring element
					// for DR boundaries
					//m_elements[e->second].neighbors[s] = 0; // 0: at least a valid value for periodic boundaries
					m_elements[e->second].neighbors[s] = m_elements.size();
				m_elements[e->second].neighborRanks[s] = m_rank;
				m_elements[e->second].boundaries[s] = boundaryCondition;

				m_mesh >> std::ws; // Skip rest of the line
			} else {
				// Not your element, ignore the boundary information
				getline(m_mesh, line);
			}
		}

		m_mesh >> std::ws;
		getline(m_mesh, line);

		if (line.find(ENDSECTION) != 0)
			abort();
	}

	void parseLocalGroups()
	{
		// This will fail if a group section come before the element section
		std::string text; // unimportant text
		m_mesh >> text;
		int groupId;
		m_mesh >> groupId;
		m_mesh >> text;
		int groupSize;
		m_mesh >> groupSize;

		std::string line;
		getline(m_mesh, line); // Skip rest of the line
		getline(m_mesh, line); // Skip group name
		getline(m_mesh, line); // Skip whatever ...

		for (int i = 0; i < groupSize; i++) {
			int element;
			m_mesh >> element;
			if (m_g2lElements.find(element-1) != m_g2lElements.end())
				m_elements[m_g2lElements[element-1]].material = groupId;
		}

		m_mesh >> std::ws;
		getline(m_mesh, line);

		if (line.find(ENDSECTION) != 0)
			abort();
	}

	void parseMPINeighborElements()
	{
#ifdef PARALLEL
		m_partition.clear();
		m_partition.seekg(0);
#endif // PARALLEL
		m_mesh.clear();
		m_mesh.seekg(m_seekElements);

		for (int i = 0; i < m_nGlobElements; i++) {
			Element element;
			element.localId = m_elements.size();
			element.rank = nextRank();

			if (element.rank != m_rank) {
				int n, t;

				m_mesh >> n; // Element number
				m_mesh >> t; // Type
				if (t != 6)
					abort();
				m_mesh >> n; // Number of points

				int commonVertices = 0;
				for (int j = 0; j < 4; j++) {
					m_mesh >> element.vertices[j];
					element.vertices[j]--;

					if (m_g2lVertices.find(element.vertices[j]) != m_g2lVertices.end())
						commonVertices++;

					// Set neighbor rank to -1, otherwise findAndUpdateNeighbors will skip the side
					element.neighbors[j] = -1;
				}

				if (commonVertices >= 3) {
					// 3 or more vertices are in our domain
					// This must be a neighbor element
					findAndUpdateNeighbors(element);

					for (int j = 0; j < 4; j++) {
						if (element.neighbors[j] >= 0) {
							// Add MPI neighbor element
							MPINeighborElement neighbor = {element.neighbors[j], element.neighborSides[j], i, j};
							m_MPINeighbors[element.rank].elements.push_back(neighbor);
						}
					}
				}

				m_mesh >> std::ws; // Skip rest of the line
			} else {
				// Ignore our elements
				std::string line;
				getline(m_mesh, line);
			}
		}

		std::string line;
		m_mesh >> std::ws;
		getline(m_mesh, line);

		if (line.find(ENDSECTION) != 0)
			abort();
	}

	void parseLocalNeighborElements()
	{
		for (unsigned int i = 0; i < m_elements.size(); i++)
			findAndUpdateNeighbors(m_elements[i]);
	}

	void parseLocalCoordinates()
	{
		m_mesh.clear();
		m_mesh.seekg(m_seekVertices);

		std::map<int, int>::iterator iter = m_g2lVertices.begin();
		for (int i = 0; i < m_nGlobVertices; i++) {
			if (iter != m_g2lVertices.end() && i == iter->first) {
				int n;

				m_mesh >> n;	// Vertex number
				for (int j = 0; j < 3; j++)
					m_mesh >> m_vertices[iter->second].coords[j];

				iter++;

				m_mesh >> std::ws; // Skip rest of the line
			} else {
				// Ignore the vertex, not in one of our elements
				std::string line;
				getline(m_mesh, line);
			}
		}

		std::string line;
		m_mesh >> std::ws;
		getline(m_mesh, line);

		if (line.find(ENDSECTION) != 0)
			abort();
	}

	void translateG2LVertices()
	{
		for (std::vector<Element>::iterator i = m_elements.begin();
				i != m_elements.end(); i++) {
			for (int j = 0; j < 4; j++)
				i->vertices[j] = m_g2lVertices[i->vertices[j]];
		}
	}

	void sortMPINeighborElements()
	{
		std::map<int, MPINeighbor>::iterator iter = m_MPINeighbors.begin();
		for (unsigned int i = 0; i < m_MPINeighbors.size(); i++) {
			iter->second.localID = i;

			if (iter->first > m_rank)
				std::sort(iter->second.elements.begin(), iter->second.elements.end(), compareLocalMPINeighbor);
			else
				std::sort(iter->second.elements.begin(), iter->second.elements.end(), compareRemoteMPINeighbor);

			// Set the MPI number of all elements
			for (unsigned int j = 0; j < iter->second.elements.size(); j++) {
				m_elements[iter->second.elements[j].localElement].mpiIndices[iter->second.elements[j].localSide]
				    = j;
			}

			iter++;
		}
	}

	/**
	 *
	 */
	void findAndUpdateNeighbors(Element &element)
	{
		// Find all elements that share one vertices with element
		std::map<int, int>::iterator potNeighbors[4];
		for (int i = 0; i < 4; i++) {
				potNeighbors[i] = m_g2lVertices.find(element.vertices[i]);
		}

		if (!(potNeighbors[0] == m_g2lVertices.end()
			|| potNeighbors[1] == m_g2lVertices.end()
			|| potNeighbors[2] == m_g2lVertices.end()
			|| element.neighbors[0] >= 0)) {
			// Find face 0 neighbor

			std::vector<int>& n0 = m_vertices[potNeighbors[0]->second].elements;
			std::vector<int>& n1 = m_vertices[potNeighbors[1]->second].elements;
			std::vector<int>& n2 = m_vertices[potNeighbors[2]->second].elements;

			std::vector<int>::const_iterator neighbor = n0.begin();
			intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end() && &element == &m_elements[*neighbor])
				// Found same element -> search for next
				intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end()) {
				updateNeighbor(element, m_elements[*neighbor], 0);
			}
		}

		if (!(potNeighbors[0] == m_g2lVertices.end()
			|| potNeighbors[1] == m_g2lVertices.end()
			|| potNeighbors[3] == m_g2lVertices.end()
			|| element.neighbors[1] >= 0)) {
			// Find face 1 neighbor

			std::vector<int>& n0 = m_vertices[potNeighbors[0]->second].elements;
			std::vector<int>& n1 = m_vertices[potNeighbors[1]->second].elements;
			std::vector<int>& n2 = m_vertices[potNeighbors[3]->second].elements;

			std::vector<int>::const_iterator neighbor = n0.begin();
			intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end() && &element == &m_elements[*neighbor])
				// Found same element -> search for next
				intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end()) {
				updateNeighbor(element, m_elements[*neighbor], 1);
			}
		}

		if (!(potNeighbors[0] == m_g2lVertices.end()
			|| potNeighbors[2] == m_g2lVertices.end()
			|| potNeighbors[3] == m_g2lVertices.end()
			|| element.neighbors[2] >= 0)) {
			// Find face 2 neighbor

			std::vector<int>& n0 = m_vertices[potNeighbors[0]->second].elements;
			std::vector<int>& n1 = m_vertices[potNeighbors[2]->second].elements;
			std::vector<int>& n2 = m_vertices[potNeighbors[3]->second].elements;

			std::vector<int>::const_iterator neighbor = n0.begin();
			intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end() && &element == &m_elements[*neighbor])
				// Found same element -> search for next
				intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end()) {
				updateNeighbor(element, m_elements[*neighbor], 2);
			}
		}

		if (!(potNeighbors[1] == m_g2lVertices.end()
			|| potNeighbors[2] == m_g2lVertices.end()
			|| potNeighbors[3] == m_g2lVertices.end()
			|| element.neighbors[3] >= 0)) {
			// Find face 3 neighbor

			std::vector<int>& n0 = m_vertices[potNeighbors[1]->second].elements;
			std::vector<int>& n1 = m_vertices[potNeighbors[2]->second].elements;
			std::vector<int>& n2 = m_vertices[potNeighbors[3]->second].elements;

			std::vector<int>::const_iterator neighbor = n0.begin();
			intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end() && &element == &m_elements[*neighbor])
				// Found same element -> search for next
				intersection(neighbor, n0.end(), n1, n2);

			if (neighbor != n0.end()) {
				updateNeighbor(element, m_elements[*neighbor], 3);
			}
		}
	}

	void updateNeighbor(Element &elem1, Element &elem2, int side1)
	{
		static const int faces[4][3] = {
				{0,2,1},
				{0,1,3},
				{0,3,2},
				{1,2,3}
		};

		// Calculate connected side of element 2
		// We use the fact that the sum of the vertices indices = side index + 3
		int side2 = -3;
		for (int i = 0; i < 3; i++) {
			side2 += std::find(elem2.vertices, elem2.vertices+4, elem1.vertices[faces[side1][i]]) - elem2.vertices;
		}

		elem1.neighbors[side1] = elem2.localId;
		elem1.neighborSides[side1] = side2;
		elem1.sideOrientations[side1] = std::find(faces[side2], faces[side2]+3,
				std::find(elem2.vertices, elem2.vertices+4, elem1.vertices[faces[side1][0]]) - elem2.vertices)
				- faces[side2];
		elem1.neighborRanks[side1] = elem2.rank;

		elem2.neighbors[side2] = elem1.localId;
		elem2.neighborSides[side2] = side1;
		elem2.sideOrientations[side2] = std::find(faces[side1], faces[side1]+3,
				std::find(elem1.vertices, elem1.vertices+4, elem2.vertices[faces[side2][0]]) - elem1.vertices)
				- faces[side1];
		elem2.neighborRanks[side2] = elem1.rank;
	}

	/**
	 * @return The next rank from the partition file or 0 if MPI not compiled with MPI
	 */
	int nextRank()
	{
		int rank;
#ifdef PARALLEL
		m_partition >> rank;
#else // PARALLEL
		rank = 0;
#endif // PARALLEL
		return rank;
	}

private:
	static const char* GAMBIT_FILE_ID;
	static const char* ENDSECTION;
	static const char* NODAL_COORDINATES;
	static const char* ELEMENT_CELLS;
	static const char* ELEMENT_GROUP;
	static const char* BOUNDARY_CONDITIONS;
};

#endif // GAMBIT_READER_H

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
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
 */

#include <cassert>
#include <cstring>

#include "SeisSol.h"
#include "WaveFieldWriter.h"
#include "Geometry/MeshReader.h"
#include "Geometry/refinement/MeshRefiner.h"

void seissol::writer::WaveFieldWriter::init(unsigned int numVars,
		int order, int numAlignedDOF,
		const MeshReader &meshReader,
		const double* dofs,  const double* pstrain, const double* integrals,
		unsigned int* map,
		int refinement, int timestep, int* outputMask, double* outputRegionBounds,
		double timeTolerance)
{
	if (!m_enabled)
		return;

	// Initialize the asynchronous module
	async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>::init();

	const int rank = seissol::MPI::mpi.rank();

	logInfo(rank) << "Initializing XDMF wave field output.";

	/** All initialization parameters */
	WaveFieldInitParam param;

	param.timestep = timestep;

	/** List of all buffer ids */
	param.bufferIds[OUTPUT_PREFIX] = addSyncBuffer(m_outputPrefix.c_str(), m_outputPrefix.size()+1, true);

	//
	// High order I/O
	//
	m_numVariables = numVars;
	m_outputFlags = new bool[numVars];
	for (size_t i = 0; i < numVars; i++)
		m_outputFlags[i] = (outputMask[i] != 0);
	// WARNING: The m_outputFlags memory might be directly used by the executor.
	// Do not modify this array after the following line
	param.bufferIds[OUTPUT_FLAGS] = addSyncBuffer(m_outputFlags, numVars*sizeof(bool), true);

	// Setup the tetrahedron refinement strategy
	refinement::TetrahedronRefiner<double>* tetRefiner = 0L;
	switch (refinement) {
	case 0:
		logInfo(rank) << "Refinement is turned off.";
		tetRefiner = new refinement::IdentityRefiner<double>();
		break;
	case 1:
		logInfo(rank) << "Refinement Startegy is \"Divide by 4\"";
		tetRefiner = new refinement::DivideTetrahedronBy4<double>();
		break;
	case 2:
		logInfo(rank) << "Refinement Startegy is \"Divide by 8\"";
		tetRefiner = new refinement::DivideTetrahedronBy8<double>();
		break;
	case 3:
		logInfo(rank) << "Refinement Startegy is \"Divide by 32\"";
		tetRefiner = new refinement::DivideTetrahedronBy32<double>();
		break;
	default:
		logError() << "Refinement Strategy is invalid!" << std::endl
			<< "Current value : " << refinement << std::endl
			<< "Valid options are :" << std::endl << "0 - No Refinement"
			<< std::endl << "1 - Refinement by 4" << std::endl
			<< "2 - Refinement by 8" << std::endl
			<< "3 - Refinement by 32";
	}

	unsigned int numElems = meshReader.getElements().size();
	unsigned int numVerts = meshReader.getVertices().size();

	// These variables are only filled if a region is to be extracted
	// Elements of the extracted region
	std::vector<const Element*> subElements;
	// The oldToNewVertexMap defines a map between old vertex index to
	// new vertex index. This is used to assign the vertex subset as well as
	// used in MeshRefiner since the elements would hold old index of vertices
	std::map<int, int> oldToNewVertexMap;
	// Vertices of the extracted region
	std::vector<const Vertex*> subVertices;
	// Mesh refiner
	refinement::MeshRefiner<double>* meshRefiner = 0L;

	// If at least one of the outputRegionBounds is non-zero then extract.
	// Otherwise use the entire region.
	// m_extractRegion = true  : Extract region
	// m_extractRegion = false : Entire region
	m_extractRegion = outputRegionBounds[0] != 0.0 ||
		outputRegionBounds[1] != 0.0 || outputRegionBounds[2] != 0.0 ||
		outputRegionBounds[3] != 0.0 || outputRegionBounds[4] != 0.0 ||
		outputRegionBounds[5] != 0.0;

	if (m_extractRegion) {
		/** Extract the elements and vertices based on user given bounds */
		// Reference to the vector containing all the elements
		const std::vector<Element>& allElements = meshReader.getElements();

		// Reference to the vector containing all the vertices
		const std::vector<Vertex>& allVertices = meshReader.getVertices();

		// m_map will store a new map from new cell index to dof index
		// the old map is contained in the "map" variable - which is a map from old
		//    cell index to dof index
		m_map = new unsigned int[numElems];

		// Extract elements based on the region specified
		for (size_t i = 0; i < numElems; i++) {
			// Store the current number of elements to check if new was added
			if (vertexInBox(outputRegionBounds, allVertices[allElements[i].vertices[0]].coords) ||
				vertexInBox(outputRegionBounds, allVertices[allElements[i].vertices[1]].coords) ||
				vertexInBox(outputRegionBounds, allVertices[allElements[i].vertices[2]].coords) ||
				vertexInBox(outputRegionBounds, allVertices[allElements[i].vertices[3]].coords)) {

				// Assign the new map
				m_map[subElements.size()] = map[i];

				// Push the address of the element into the vector
				subElements.push_back(&(allElements[i]));

				// Push the vertices into the map which makes sure that the entries are unique
				oldToNewVertexMap.insert(std::pair<int,int>(allElements[i].vertices[0], oldToNewVertexMap.size()));
				oldToNewVertexMap.insert(std::pair<int,int>(allElements[i].vertices[1], oldToNewVertexMap.size()));
				oldToNewVertexMap.insert(std::pair<int,int>(allElements[i].vertices[2], oldToNewVertexMap.size()));
				oldToNewVertexMap.insert(std::pair<int,int>(allElements[i].vertices[3], oldToNewVertexMap.size()));
			}
		}

		subVertices.resize(oldToNewVertexMap.size());

		// Loop over the map and assign the vertices
		for (std::map<int,int>::iterator it=oldToNewVertexMap.begin(); it!=oldToNewVertexMap.end(); ++it)
			subVertices.at(it->second) = &allVertices[it->first];

		numElems = subElements.size();
		numVerts = subVertices.size();

		meshRefiner = new refinement::MeshRefiner<double>(subElements, subVertices,
			oldToNewVertexMap, *tetRefiner);
	} else {
		meshRefiner = new refinement::MeshRefiner<double>(meshReader, *tetRefiner);
	}

	logInfo(rank) << "Refinement class initialized";
	logDebug() << "Cells : "
			<< numElems << "refined-to ->"
			<< meshRefiner->getNumCells();
	logDebug() << "Vertices : "
			<< numVerts << "refined-to ->"
			<< meshRefiner->getNumVertices();

	// Initialize the variable subsampler
	m_variableSubsampler = new refinement::VariableSubsampler<double>(
			numElems, *tetRefiner, order, numVars, numAlignedDOF);

	logInfo(rank) << "VariableSubsampler initialized";

	// Delete the tetRefiner since it is no longer required
	delete tetRefiner;

	// Cells are a bit complicated because the vertex filter will now longer work if we just use the buffer
	// We will add the offset later
	const unsigned int* const_cells;
#ifdef USE_MPI
	// Add the offset to the cells
	MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().groupComm();
	unsigned int offset = meshRefiner->getNumVertices();
	MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
	offset -= meshRefiner->getNumVertices();

	// Add the offset to all cells
	unsigned int* cells = new unsigned int[meshRefiner->getNumCells() * 4];
	for (unsigned int i = 0; i < meshRefiner->getNumCells() * 4; i++)
		cells[i] = meshRefiner->getCellData()[i] + offset;
	const_cells = cells;
#else // USE_MPI
	const_cells = meshRefiner->getCellData();
#endif // USE_MPI

	// Create mesh buffers
	param.bufferIds[CELLS] = addSyncBuffer(const_cells, meshRefiner->getNumCells() * 4 * sizeof(unsigned int));
	param.bufferIds[VERTICES] = addSyncBuffer(meshRefiner->getVertexData(), meshRefiner->getNumVertices() * 3 * sizeof(double));

	// Create data buffers
	bool first = false;
	for (unsigned int i = 0; i < numVars; i++) {
		if (m_outputFlags[i]) {
			unsigned int id = addBuffer(0L, meshRefiner->getNumCells() * sizeof(double));
			if (!first) {
				param.bufferIds[VARIABLE0] = id;
				first = true;
			}
		}
	}

	// Save number of cells
	m_numCells = meshRefiner->getNumCells();

	// Set up for low order output flags
	m_lowOutputFlags = new bool[WaveFieldWriterExecutor::NUM_TOTALLOWVARS];
	for (size_t i = 0; i < WaveFieldWriterExecutor::NUM_LOWVARIABLES; i++) {
		m_lowOutputFlags[i] = (pstrain != 0L);
	}
	for (size_t i = 0; i < WaveFieldWriterExecutor::NUM_INTEGRATED_VARIABLES; i++) {
		m_lowOutputFlags[i+WaveFieldWriterExecutor::NUM_LOWVARIABLES] = seissol::SeisSol::main.postProcessor().getIntegrationMask()[i];
	}
	param.bufferIds[LOW_OUTPUT_FLAGS] = addSyncBuffer(m_lowOutputFlags, WaveFieldWriterExecutor::NUM_TOTALLOWVARS*sizeof(bool), true);
	m_numIntegratedVariables = seissol::SeisSol::main.postProcessor().getNumberOfVariables();

	//
	//  Low order I/O
	//
	refinement::MeshRefiner<double>* pLowMeshRefiner = 0L;
#ifdef USE_MPI
	unsigned int* lowCells = 0L;
#endif // USE_MPI
	const unsigned int* const_lowCells;
	if (pstrain || integrals) {
		logInfo(rank) << "Initialize low order output";

		// Refinement strategy (no refinement)
		refinement::IdentityRefiner<double> lowTetRefiner;

		// Mesh refiner
		if (m_extractRegion) {
			pLowMeshRefiner = new refinement::MeshRefiner<double>(subElements, subVertices,
				oldToNewVertexMap, lowTetRefiner);
		} else {
			pLowMeshRefiner = new refinement::MeshRefiner<double>(meshReader, lowTetRefiner);
		}

#ifdef USE_MPI
		// Same offset issue for the normal mesh
		MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().groupComm();
		unsigned int offset = meshRefiner->getNumVertices();
		MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
		offset -= meshRefiner->getNumVertices();

		// Add the offset to all cells
		lowCells = new unsigned int[meshRefiner->getNumCells() * 4];
		for (unsigned int i = 0; i < meshRefiner->getNumCells() * 4; i++)
			cells[i] = meshRefiner->getCellData()[i] + offset;
		const_lowCells = cells;
#else // USE_MPI
		const_lowCells = meshRefiner->getCellData();
#endif // USE_MPI

		// Create mesh buffers
		param.bufferIds[LOWCELLS] = addSyncBuffer(const_lowCells,
			pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
		param.bufferIds[LOWVERTICES] = addSyncBuffer(pLowMeshRefiner->getVertexData(),
			pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));

		// Create data buffers
		param.bufferIds[LOWVARIABLE0] = addBuffer(0L, pLowMeshRefiner->getNumCells() * sizeof(double));
		int numLowVars = 0;
		if (pstrain && !integrals) {
			numLowVars = WaveFieldWriterExecutor::NUM_LOWVARIABLES;
		} else if (integrals && !pstrain) {
			numLowVars = m_numIntegratedVariables;
		} else {
			numLowVars = WaveFieldWriterExecutor::NUM_LOWVARIABLES + m_numIntegratedVariables;
		}
		logInfo(rank) << "Total number of extra variables " << numLowVars;

		for (unsigned int i = 1; i < numLowVars; i++)
			addBuffer(0L, pLowMeshRefiner->getNumCells() * sizeof(double));

		// Save number of cells
		m_numLowCells = pLowMeshRefiner->getNumCells();
	} else {
		// No low order output
		param.bufferIds[LOWCELLS] = -1;
		param.bufferIds[LOWVERTICES] = -1;
		param.bufferIds[LOWVARIABLE0] = -1;
	}

	//
	// Send all buffers for initialization
	//
	sendBuffer(param.bufferIds[OUTPUT_PREFIX], m_outputPrefix.size()+1);

	sendBuffer(param.bufferIds[OUTPUT_FLAGS], numVars*sizeof(bool));

	sendBuffer(param.bufferIds[CELLS], meshRefiner->getNumCells() * 4 * sizeof(unsigned int));
	sendBuffer(param.bufferIds[VERTICES], meshRefiner->getNumVertices() * 3 * sizeof(double));

	if (pstrain || integrals) {
		sendBuffer(param.bufferIds[LOWCELLS], pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
		sendBuffer(param.bufferIds[LOWVERTICES], pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));
		sendBuffer(param.bufferIds[LOW_OUTPUT_FLAGS], WaveFieldWriterExecutor::NUM_TOTALLOWVARS*sizeof(bool));
	}

	// Initialize the executor
	callInit(param);

	// Remove buffers
	removeBuffer(param.bufferIds[OUTPUT_PREFIX]);
	removeBuffer(param.bufferIds[CELLS]);
	removeBuffer(param.bufferIds[VERTICES]);
	if (pstrain || integrals) {
		removeBuffer(param.bufferIds[LOWCELLS]);
		removeBuffer(param.bufferIds[LOWVERTICES]);
	}

	// Remove the low mesh refiner if it was setup
	if (pLowMeshRefiner)
		delete pLowMeshRefiner;

#ifdef USE_MPI
	delete [] cells;
	delete [] lowCells;
#endif // USE_MPI

	// Save dof/map pointer
	m_dofs = dofs;
	m_pstrain = pstrain;
	m_integrals = integrals;
	if (!m_extractRegion) {
		m_map = map;
	}

	m_timestep = timestep;
	m_variableBufferIds[0] = param.bufferIds[VARIABLE0];
	m_variableBufferIds[1] = param.bufferIds[LOWVARIABLE0];

	delete meshRefiner;
}

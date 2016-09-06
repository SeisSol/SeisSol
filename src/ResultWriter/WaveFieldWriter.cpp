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
		const double* dofs,  const double* pstrain,
		const unsigned int* map,
		int refinement, int timestep, int* outputMask,
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

	// Refine the mesh
	refinement::MeshRefiner<double> meshRefiner(meshReader, *tetRefiner);

	logInfo(rank) << "Refinement class initialized";
	logDebug() << "Cells : "
			<< meshReader.getElements().size() << "refined-to ->"
			<< meshRefiner.getNumCells();
	logDebug() << "Vertices : "
			<< meshReader.getVertices().size()
			<< "refined-to ->"
			<< meshRefiner.getNumVertices();

	// Initialize the variable subsampler
	m_variableSubsampler = new refinement::VariableSubsampler<double>(
			meshReader.getElements().size(),
			*tetRefiner, order, numVars, numAlignedDOF);

	logInfo(rank) << "VariableSubsampler initialized";

	// Delete the tetRefiner since it is no longer required
	delete tetRefiner;

	// Cells are a bit complicated because the vertex filter will now longer work if we just use the buffer
	// We will add the offset later
	const unsigned int* const_cells;
#ifdef USE_MPI
	// Add the offset to the cells
	MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().groupComm();
	unsigned int offset = meshRefiner.getNumVertices();
	MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
	offset -= meshRefiner.getNumVertices();

	// Add the offset to all cells
	unsigned int* cells = new unsigned int[meshRefiner.getNumCells() * 4];
	for (unsigned int i = 0; i < meshRefiner.getNumCells() * 4; i++)
		cells[i] = meshRefiner.getCellData()[i] + offset;
	const_cells = cells;
#else // USE_MPI
	const_cells = meshRefiner.getCellData();
#endif // USE_MPI

	// Create mesh buffers
	param.bufferIds[CELLS] = addSyncBuffer(const_cells, meshRefiner.getNumCells() * 4 * sizeof(unsigned int));
	param.bufferIds[VERTICES] = addSyncBuffer(meshRefiner.getVertexData(), meshRefiner.getNumVertices() * 3 * sizeof(double));

	// Create data buffers
	bool first = false;
	for (unsigned int i = 0; i < numVars; i++) {
		if (m_outputFlags[i]) {
			unsigned int id = addBuffer(0L, meshRefiner.getNumCells() * sizeof(double));
			if (!first) {
				param.bufferIds[VARIABLE0] = id;
				first = true;
			}
		}
	}

	// Save number of cells
	m_numCells = meshRefiner.getNumCells();

	//
	//  Low order I/O
	//
	refinement::MeshRefiner<double>* pLowMeshRefiner = 0L;
#ifdef USE_MPI
	unsigned int* lowCells = 0L;
#endif // USE_MPI
	const unsigned int* const_lowCells;
	if (pstrain) {
		logInfo(rank) << "Initialize low order output";

		// Refinement strategy (no refinement)
		refinement::IdentityRefiner<double> lowTetRefiner;

		// Mesh refiner
		pLowMeshRefiner = new refinement::MeshRefiner<double>(meshReader, lowTetRefiner);

#ifdef USE_MPI
		// Same offset issue for the normal mesh
		MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().groupComm();
		unsigned int offset = meshRefiner.getNumVertices();
		MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
		offset -= meshRefiner.getNumVertices();

		// Add the offset to all cells
		lowCells = new unsigned int[meshRefiner.getNumCells() * 4];
		for (unsigned int i = 0; i < meshRefiner.getNumCells() * 4; i++)
			cells[i] = meshRefiner.getCellData()[i] + offset;
		const_lowCells = cells;
#else // USE_MPI
		const_lowCells = meshRefiner.getCellData();
#endif // USE_MPI

		// Create mesh buffers
		param.bufferIds[LOWCELLS] = addSyncBuffer(const_lowCells,
			pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
		param.bufferIds[LOWVERTICES] = addSyncBuffer(pLowMeshRefiner->getVertexData(),
			pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));

		// Create data buffers
		param.bufferIds[LOWVARIABLE0] = addBuffer(0L, pLowMeshRefiner->getNumCells() * sizeof(double));
		for (unsigned int i = 1; i < WaveFieldWriterExecutor::NUM_LOWVARIABLES; i++)
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

	sendBuffer(param.bufferIds[CELLS], meshRefiner.getNumCells() * 4 * sizeof(unsigned int));
	sendBuffer(param.bufferIds[VERTICES], meshRefiner.getNumVertices() * 3 * sizeof(double));

	if (pstrain) {
		sendBuffer(param.bufferIds[LOWCELLS], pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
		sendBuffer(param.bufferIds[LOWVERTICES], pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));
	}

	// Initialize the executor
	callInit(param);

	// Remove buffers
	removeBuffer(param.bufferIds[OUTPUT_PREFIX]);
	removeBuffer(param.bufferIds[CELLS]);
	removeBuffer(param.bufferIds[VERTICES]);
	if (pstrain) {
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
	m_map = map;

	m_timestep = timestep;
	m_variableBufferIds[0] = param.bufferIds[VARIABLE0];
	m_variableBufferIds[1] = param.bufferIds[LOWVARIABLE0];
}

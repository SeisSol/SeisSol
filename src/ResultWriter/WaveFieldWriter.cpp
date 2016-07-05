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
		int refinement, int timestep,
		double timeTolerance)
{
	const int rank = seissol::MPI::mpi.rank();

	// Set the executor
	setUp();

	if (!m_enabled)
		return;

	logInfo(rank) << "Initializing XDMF wave field output.";

	if (m_waveFieldWriter != 0L)
		logError() << "Wave field writer already initialized";

	// Buffer for the mesh file name
#ifdef USE_ASYNC_MPI
	m_bufferIds[OUTPUT_PREFIX] = addBuffer(m_outputPrefix.size()+1);
#endif // USE_ASYNC_MPI

	//
	// High order I/O
	//
	if (static_cast<unsigned int>(numVars) != MAX_VARIABLES)
		logError()
				<< "XDMF output supports currently only" << MAX_VARIABLES
				<< "variables. Number of variables specified:" << numVars;
	m_numVariables = numVars;

	// Currently all variables have to be chosen.
	m_outputFlags.resize(numVars);
	std::fill(m_outputFlags.begin(), m_outputFlags.end(), true);

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

	m_numCells = meshRefiner.getNumCells();

	// Create the buffers (for mesh and data)
#ifdef USE_ASYNC_MPI
	m_bufferIds[CELLS] = addBuffer(m_numCells * 4 * sizeof(unsigned int));
	m_bufferIds[VERTICES] = addBuffer(meshRefiner.getNumVertices() * 3 * sizeof(double));
#endif // USE_ASYNC_MPI

#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
	for (int i = 0; i < MAX_VARIABLES; i++) {
		if (m_outputFlags[i])
			m_bufferIds[VARIABLE0+i] = addBuffer(m_numCells * sizeof(double));
	}
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)

	// Local buffer to extract the data
	m_outputBuffer = new double[meshRefiner.getNumCells()];
	assert(meshRefiner.getNumCells() >= meshReader.getElements().size());

	//
	//  Low order I/O
	//
	refinement::MeshRefiner<double>* pLowMeshRefiner = 0L;
	if (pstrain) {
		logInfo(rank) << "Initialize low order output";

		// Refinement strategy (no refinement)
		refinement::IdentityRefiner<double> lowTetRefiner;

		// Mesh refiner
		pLowMeshRefiner = new refinement::MeshRefiner<double>(meshReader, lowTetRefiner);

		m_numLowCells = pLowMeshRefiner->getNumCells();

		// Create the buffer (for mesh and data)
#ifdef USE_ASYNC_MPI
		m_bufferIds[LOWCELLS] = addBuffer(m_numLowCells * 4 * sizeof(unsigned int));
		m_bufferIds[LOWVERTICES] = addBuffer(pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));
#endif // USE_ASYNC_MPI

#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
			for (int i = 0; i < MAX_LOWVARIABLES; i++)
				m_bufferIds[LOWVARIABLE0+i] = addBuffer(m_numLowCells * sizeof(double));
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
		}

	// Save dof/map pointer
	m_dofs = dofs;
	m_pstrain = pstrain;
	m_map = map;

	m_timestep = timestep;

#ifdef USE_ASYNC_MPI
	// Fill the mesh buffers
	fillBuffer(OUTPUT_PREFIX, m_outputPrefix.c_str(), m_outputPrefix.size()+1);

	// Cells are a bit complicated because the vertex filter will now longer work if we just use the buffer
	MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().scheduler().groupComm();
	unsigned int offset = meshRefiner.getNumVertices();
	MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
	offset -= meshRefiner.getNumVertices();

	// Add the offset to all cells
	unsigned int* cells = new unsigned int[meshRefiner.getNumCells() * 4];
	for (unsigned int i = 0; i < meshRefiner.getNumCells() * 4; i++)
		cells[i] = meshRefiner.getCellData()[i] + offset;

	fillBuffer(m_bufferIds[CELLS], cells, meshRefiner.getNumCells() * 4 * sizeof(unsigned int));
	fillBuffer(m_bufferIds[VERTICES], meshRefiner.getVertexData(), meshRefiner.getNumVertices() * 3 * sizeof(double));
	if (pstrain) {
		// Same offset issue for the low order mesh
		offset = pLowMeshRefiner->getNumVertices();
		MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
		offset -= pLowMeshRefiner->getNumVertices();

		for (unsigned int i = 0; i < pLowMeshRefiner->getNumCells() * 4; i++)
			cells[i] = pLowMeshRefiner->getCellData()[i] + offset;

		fillBuffer(m_bufferIds[LOWCELLS], cells,
				pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
		fillBuffer(m_bufferIds[LOWVERTICES], pLowMeshRefiner->getVertexData(),
				pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));
	}

	delete [] cells;
#endif // USE_ASYNC_MPI

	// Initialize the executor
	WaveFieldInitParam param;
	param.timestep = timestep;
	param.numVars = numVars;
#ifdef USE_ASYNC_MPI
	assert(sizeof(param.bufferIds) == sizeof(m_bufferIds));
	memcpy(param.bufferIds, m_bufferIds, sizeof(m_bufferIds));
#else // USE_ASYNC_MPI
	param.numCells = meshRefiner.getNumCells();
	param.cells = meshRefiner.getCellData();
	param.numVertices = meshRefiner.getNumVertices();
	param.vertices = meshRefiner.getVertexData();

	if (pstrain) {
		param.numLowCells = pLowMeshRefiner->getNumCells();
		param.lowCells = pLowMeshRefiner->getCellData();
		param.numLowVertices = pLowMeshRefiner->getNumVertices();
		param.lowVertices = pLowMeshRefiner->getVertexData();
	} else {
		param.numLowCells = 0;
		param.lowCells = 0L;
		param.numLowVertices = 0;
		param.lowVertices = 0L;
	}
#endif // USE_ASYNC_MPI

	callInit(param);

	// Remove the low mesh refiner if it was setup
	if (pLowMeshRefiner)
		delete pLowMeshRefiner;
}


void seissol::writer::WaveFieldWriter::execInit(const WaveFieldInitParam &param)
{
	// Enable output on the executor
	enable();

	const int rank = seissol::MPI::mpi.rank();

	//
	// High order I/O
	//
	m_numVariables = param.numVars;
	assert(m_numVariables == 9); // Currently nothing else is supported
	std::vector<const char*> variables(m_numVariables);
	variables[0] = "sigma_xx";
	variables[1] = "sigma_yy";
	variables[2] = "sigma_zz";
	variables[3] = "sigma_xy";
	variables[4] = "sigma_yz";
	variables[5] = "sigma_xz";
	variables[6] = "u";
	variables[7] = "v";
	variables[8] = "w";

	// Initialize the I/O handler and write the mesh
	const char* outputPrefix = m_outputPrefix.c_str();
#ifdef USE_ASYNC_MPI
	memcpy(m_bufferIds, param.bufferIds, sizeof(m_bufferIds));
	outputPrefix = static_cast<const char*>(buffer(m_bufferIds[OUTPUT_PREFIX]));
#endif // USE_ASYNC_MPI
	m_waveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
			rank, outputPrefix, variables, param.timestep);

#ifdef USE_MPI
	MPI_Comm_dup(seissol::MPI::mpi.comm(), &m_comm);
	m_waveFieldWriter->setComm(m_comm);
#endif // USE_MPI

#ifdef USE_ASYNC_MPI
	size_t numCells = bufferSize(m_bufferIds[CELLS]) / (4*sizeof(unsigned int));
	const unsigned int* cells = static_cast<const unsigned int*>(buffer(m_bufferIds[CELLS]));
	size_t numVertices = bufferSize(m_bufferIds[VERTICES]) / (3*sizeof(double));
	const double* vertices = static_cast<const double*>(buffer(m_bufferIds[VERTICES]));
#else // USE_ASYNC_MPI
	size_t numCells = param.numCells;
	const unsigned int* cells = param.cells;
	size_t numVertices = param.numVertices;
	const double* vertices = param.vertices;
#endif // USE_ASYNC_MPI
	m_waveFieldWriter->init(numCells, cells, numVertices, vertices, true);

	logInfo(rank) << "High order output initialized";

	//
	// Low order I/O
	//
#ifdef USE_ASYNC_MPI
	if (m_bufferIds[LOWCELLS] >= 0) {
		numCells = bufferSize(m_bufferIds[LOWCELLS]) / (4*sizeof(unsigned int));
		cells = static_cast<const unsigned int*>(buffer(m_bufferIds[LOWCELLS]));
		numVertices = bufferSize(m_bufferIds[LOWVERTICES]) / (3*sizeof(double));
		vertices = static_cast<const double*>(buffer(m_bufferIds[LOWVERTICES]));
	} else {
		numCells = 0; // Low order output not enabled
	}
#else // USE_ASYNC_MPI
	numCells = param.numLowCells;
	cells = param.lowCells;
	numVertices = param.numLowVertices;
	vertices = param.lowVertices;
#endif // USE_ASYNC_MPI

	if (numCells > 0) {
		// Pstrain enabled

		// Variables
		std::vector<const char*> lowVariables(MAX_LOWVARIABLES);
		lowVariables[0] = "ep_xx";
		lowVariables[1] = "ep_yy";
		lowVariables[2] = "ep_zz";
		lowVariables[3] = "ep_xy";
		lowVariables[4] = "ep_yz";
		lowVariables[5] = "ep_xz";
		lowVariables[6] = "eta";

		m_lowWaveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
				rank, (std::string(outputPrefix)+"-low").c_str(), lowVariables, param.timestep);

		logInfo(rank) << "Low order output initialized";

#ifdef USE_MPI
		m_lowWaveFieldWriter->setComm(m_comm);
#endif // USE_MPI

		m_lowWaveFieldWriter->init(numCells, cells, numVertices, vertices, true);
	}

	logInfo(rank) << "Initializing XDMF wave field output. Done.";
}

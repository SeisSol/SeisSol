/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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
#include "Monitoring/instrumentation.fpp"
#include <Modules/Modules.h>

void seissol::writer::WaveFieldWriter::setUp()
{
  setExecutor(m_executor);
  if (isAffinityNecessary()) {
    const auto freeCpus = SeisSol::main.getPinning().getFreeCPUsMask();
    logInfo(seissol::MPI::mpi.rank()) << "Wave field writer thread affinity:"
      << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
    setAffinityIfNecessary(freeCpus);
  }
}

void seissol::writer::WaveFieldWriter::enable()
{
	m_enabled = true;
	seissol::SeisSol::main.checkPointManager().header().add(m_timestepComp);
}

seissol::refinement::TetrahedronRefiner<double>* seissol::writer::WaveFieldWriter::createRefiner(int refinement) {
  int const rank = seissol::MPI::mpi.rank();
  refinement::TetrahedronRefiner<double>* tetRefiner = 0L;
  switch (refinement) {
	case 0:
		logInfo(rank) << "Refinement is turned off.";
		tetRefiner = new refinement::IdentityRefiner<double>();
		break;
	case 1:
		logInfo(rank) << "Refinement Strategy is \"Divide by 4\"";
		tetRefiner = new refinement::DivideTetrahedronBy4<double>();
		break;
	case 2:
		logInfo(rank) << "Refinement Strategy is \"Divide by 8\"";
		tetRefiner = new refinement::DivideTetrahedronBy8<double>();
		break;
	case 3:
		logInfo(rank) << "Refinement Strategy is \"Divide by 32\"";
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
  return tetRefiner;
}

unsigned const* seissol::writer::WaveFieldWriter::adjustOffsets(refinement::MeshRefiner<double>* meshRefiner) {
  unsigned const* const_cells;
// Cells are a bit complicated because the vertex filter will now longer work if we just use the buffer
// We will add the offset later
#ifdef USE_MPI
	// Add the offset to the cells
	MPI_Comm groupComm = seissol::SeisSol::main.asyncIO().groupComm();
	unsigned int offset = meshRefiner->getNumVertices();
	MPI_Scan(MPI_IN_PLACE, &offset, 1, MPI_UNSIGNED, MPI_SUM, groupComm);
	offset -= meshRefiner->getNumVertices();

	// Add the offset to all cells
	unsigned int* cells = new unsigned int[meshRefiner->getNumCells() * 4];
	for (unsigned int i = 0; i < meshRefiner->getNumCells() * 4; i++) {
		cells[i] = meshRefiner->getCellData()[i] + offset;
  }
	const_cells = cells;
#else // USE_MPI
	const_cells = meshRefiner->getCellData();
#endif // USE_MPI
  return const_cells;
}

void seissol::writer::WaveFieldWriter::init(unsigned int numVars,
    int order, int numAlignedDOF,
    const MeshReader &meshReader, const std::vector<unsigned> &LtsClusteringData,
    const real* dofs,  const real* pstrain, const real* integrals,
    unsigned int* map,
    int refinement, int* outputMask, double* outputRegionBounds,
    xdmfwriter::BackendType backend)
{
	if (!m_enabled)
		return;

	// Initialize the asynchronous module
	async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>::init();

  Modules::registerHook(*this, SIMULATION_START);
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);

	const int rank = seissol::MPI::mpi.rank();

	logInfo(rank) << "Initializing XDMF wave field output.";

	/** All initialization parameters */
	WaveFieldInitParam param;

	param.timestep = seissol::SeisSol::main.checkPointManager().header().value(m_timestepComp);

	/** List of all buffer ids */
	param.bufferIds[OUTPUT_PREFIX] = addSyncBuffer(m_outputPrefix.c_str(), m_outputPrefix.size()+1, true);

  param.backend = backend;

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
	refinement::TetrahedronRefiner<double>* tetRefiner = createRefiner(refinement);

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

	const unsigned int* const_cells = adjustOffsets(meshRefiner);

	// Create mesh buffers
	param.bufferIds[CELLS] = addSyncBuffer(const_cells, meshRefiner->getNumCells() * 4 * sizeof(unsigned int));
	param.bufferIds[VERTICES] = addSyncBuffer(meshRefiner->getVertexData(), meshRefiner->getNumVertices() * 3 * sizeof(double));

	// Create data buffers
	bool first = false;
	for (unsigned int i = 0; i < numVars; i++) {
		if (m_outputFlags[i]) {
			unsigned int id = addBuffer(0L, meshRefiner->getNumCells() * sizeof(real));
			if (!first) {
				param.bufferIds[VARIABLE0] = id;
				first = true;
			}
		}
	}

	// Save number of cells
	m_numCells = meshRefiner->getNumCells();
	// Set up for low order output flags
	m_lowOutputFlags = new bool[WaveFieldWriterExecutor::NUM_LOWVARIABLES];
	m_numIntegratedVariables = seissol::SeisSol::main.postProcessor().getNumberOfVariables();
	for (size_t i = 0; i < WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES; i++) {
		m_lowOutputFlags[i] = (pstrain != 0L);
	}
	seissol::SeisSol::main.postProcessor().getIntegrationMask(&m_lowOutputFlags[WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES]);
	param.bufferIds[LOW_OUTPUT_FLAGS] = addSyncBuffer(m_lowOutputFlags, WaveFieldWriterExecutor::NUM_LOWVARIABLES*sizeof(bool), true);
	//
	//  Low order I/O
	//
	refinement::MeshRefiner<double>* pLowMeshRefiner = 0L;
	const unsigned int* const_lowCells = 0L;
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

    const_lowCells = adjustOffsets(pLowMeshRefiner);

		// Create mesh buffers
		param.bufferIds[LOWCELLS] = addSyncBuffer(const_lowCells,
			pLowMeshRefiner->getNumCells() * 4 * sizeof(unsigned int));
		param.bufferIds[LOWVERTICES] = addSyncBuffer(pLowMeshRefiner->getVertexData(),
			pLowMeshRefiner->getNumVertices() * 3 * sizeof(double));

		// Create data buffers
		param.bufferIds[LOWVARIABLE0] = addBuffer(0L, pLowMeshRefiner->getNumCells() * sizeof(real));
		int numLowVars = 0;

		if (pstrain) numLowVars += WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES;

		if (integrals) numLowVars += m_numIntegratedVariables;

		for (int i = 1; i < numLowVars; i++)
			addBuffer(0L, pLowMeshRefiner->getNumCells() * sizeof(real));

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
		sendBuffer(param.bufferIds[LOW_OUTPUT_FLAGS], WaveFieldWriterExecutor::NUM_LOWVARIABLES*sizeof(bool));
	}

	// Initialize the executor
	callInit(param);
  m_executor.setClusteringData(LtsClusteringData.data());

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
	delete [] const_cells;
	delete [] const_lowCells;
#endif // USE_MPI

	// Save dof/map pointer
	m_dofs = dofs;
	m_pstrain = pstrain;
	m_integrals = integrals;
	if (!m_extractRegion) {
		m_map = map;
	}

	m_variableBufferIds[0] = param.bufferIds[VARIABLE0];
	m_variableBufferIds[1] = param.bufferIds[LOWVARIABLE0];

	delete meshRefiner;
}

void seissol::writer::WaveFieldWriter::write(double time)
{
	SCOREP_USER_REGION("WaveFieldWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	if (!m_enabled)
		return;

	m_stopwatch.start();

	const int rank = seissol::MPI::mpi.rank();

	SCOREP_USER_REGION_DEFINE(r_wait);
	SCOREP_USER_REGION_BEGIN(r_wait, "wavfieldwriter_wait", SCOREP_USER_REGION_TYPE_COMMON);
	logInfo(rank) << "Waiting for last wave field.";
	wait();
	SCOREP_USER_REGION_END(r_wait);

	logInfo(rank) << "Writing wave field at time" << utils::nospace <<  time << '.';

	unsigned int nextId = m_variableBufferIds[0];
	for (unsigned int i = 0; i < m_numVariables; i++) {
		if (!m_outputFlags[i])
			continue;

		real* managedBuffer = async::Module<WaveFieldWriterExecutor,
				WaveFieldInitParam, WaveFieldParam>::managedBuffer<real*>(nextId);
		m_variableSubsampler->get(m_dofs, m_map, i, managedBuffer);

		sendBuffer(nextId, m_numCells*sizeof(real));

		nextId++;
	}

	// nextId is required in a manner similar to above for writing integrated variables
	nextId = 0;
	if (m_pstrain) {
		for (unsigned int i = 0; i < WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES; i++) {
			real* managedBuffer = async::Module<WaveFieldWriterExecutor,
					WaveFieldInitParam, WaveFieldParam>::managedBuffer<real*>(m_variableBufferIds[1]+i);

#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
#endif // _OPENMP
			for (unsigned int j = 0; j < m_numLowCells; j++)
				managedBuffer[j] = m_pstrain[m_map[j]
						* WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES + i];

			sendBuffer(m_variableBufferIds[1]+i, m_numLowCells*sizeof(real));
		}
		nextId = WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES;
	}

	// This offset is used to access the correct variable in m_integrals
	// If pstrain is enabled then the offset is set to NUM_PLASTICITY_VARIABLES otherwise it is set to 0
	unsigned int offset = nextId;

	if (m_integrals) {
		for (unsigned int i = 0; i < WaveFieldWriterExecutor::NUM_INTEGRATED_VARIABLES; i++) {
			if (!m_lowOutputFlags[i+WaveFieldWriterExecutor::NUM_PLASTICITY_VARIABLES])
				continue;
			real* managedBuffer = async::Module<WaveFieldWriterExecutor,
			WaveFieldInitParam, WaveFieldParam>::managedBuffer<real*>(m_variableBufferIds[1]+nextId);

#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
#endif // _OPENMP
			for (unsigned int j = 0; j < m_numLowCells; j++)
				managedBuffer[j] = m_integrals[m_map[j]
						* m_numIntegratedVariables + nextId - offset];

			sendBuffer(m_variableBufferIds[1]+nextId, m_numLowCells*sizeof(real));
			nextId++;
		}
	}

	WaveFieldParam param;
	param.time = time;
	call(param);

	// Update last time step
	seissol::SeisSol::main.checkPointManager().header().value(m_timestepComp)++;

	m_stopwatch.pause();

	logInfo(rank) << "Writing wave field at time" << utils::nospace << time << ". Done.";
}

void seissol::writer::WaveFieldWriter::simulationStart()
{
	syncPoint(0.0);
}

void seissol::writer::WaveFieldWriter::syncPoint(double currentTime)
{
	write(currentTime);
}

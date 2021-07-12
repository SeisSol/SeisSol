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

#ifndef WAVE_FIELD_WRITER_EXECUTOR_H
#define WAVE_FIELD_WRITER_EXECUTOR_H

#include "Parallel/MPI.h"

#include <cassert>
#include <vector>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "async/ExecInfo.h"

#include "Monitoring/Stopwatch.h"

namespace seissol
{

namespace writer
{

/** Buffer tags for asynchronous IO */
enum BufferTags {
	OUTPUT_PREFIX,
	OUTPUT_FLAGS,
	CELLS,
	VERTICES,
	VARIABLE0,
	LOWCELLS,
	LOWVERTICES,
	LOW_OUTPUT_FLAGS,
	LOWVARIABLE0,
	BUFFERTAG_MAX = LOWVARIABLE0
};

struct WaveFieldInitParam
{
	int timestep;

	int bufferIds[BUFFERTAG_MAX+1];
  xdmfwriter::BackendType backend;
};

struct WaveFieldParam
{
	double time;
};

class WaveFieldWriterExecutor
{
private:
	/** The XMDF Writer used for the wave field */
	xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>* m_waveFieldWriter;

	/** The XDMF Writer for low order data */
	xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>* m_lowWaveFieldWriter;

	/** Buffer id for the first variable for high and low order output */
	unsigned int m_variableBufferIds[2];

	/** The total number of (high order) variables */
	unsigned int m_numVariables;

	/** Flag indicated which variables should be written */
	const bool* m_outputFlags;

	/** Flags indicating which low order variables should be written */
	const bool* m_lowOutputFlags;

#ifdef USE_MPI
	/** The MPI communicator for the XDMF writer */
	MPI_Comm m_comm;
#endif // USE_MPI

	/** Stopwatch for the wave field backend */
	Stopwatch m_stopwatch;

public:
	WaveFieldWriterExecutor()
		: m_waveFieldWriter(0L),
		  m_lowWaveFieldWriter(0L),
		  m_numVariables(0),
		  m_outputFlags(0L),
		  m_lowOutputFlags(0L)
#ifdef USE_MPI
		  , m_comm(MPI_COMM_NULL)
#endif // USE_MPI
	{ }

	/**
	 * Initialize the XDMF writers
	 */
	void execInit(const async::ExecInfo &info, const WaveFieldInitParam &param)
	{
		if (m_waveFieldWriter != 0L)
			logError() << "Wave field writer already initialized";

		int rank = seissol::MPI::mpi.rank();

		xdmfwriter::BackendType type = param.backend;

		const char* outputPrefix = static_cast<const char*>(info.buffer(param.bufferIds[OUTPUT_PREFIX]));

		//
		// High order I/O
		//
		m_numVariables = info.bufferSize(param.bufferIds[OUTPUT_FLAGS]) / sizeof(bool);
		m_outputFlags = static_cast<const bool*>(info.buffer(param.bufferIds[OUTPUT_FLAGS]));

		const char* varNames[9] = {
			"sigma_xx",
			"sigma_yy",
			"sigma_zz",
			"sigma_xy",
			"sigma_yz",
			"sigma_xz",
			"u",
			"v",
			"w",
		};

		std::vector<const char*> variables;
		for (unsigned int i = 0; i < m_numVariables; i++) {
			if (m_outputFlags[i]) {
        assert(i < 9);
				variables.push_back(varNames[i]);
      }
		}

#ifdef USE_MPI
	// Split the communicator into two - those containing vertices and those
	//  not containing any vertices.
	int commColour = (info.bufferSize(param.bufferIds[CELLS]) == 0)?0:1;
	MPI_Comm_split(seissol::MPI::mpi.comm(), commColour, rank, &m_comm);
	// Start the if statement
	if (info.bufferSize(param.bufferIds[CELLS]) != 0) {
		// Get the new rank
		MPI_Comm_rank(m_comm, &rank);
#endif // USE_MPI

		// Initialize the I/O handler and write the mesh
		m_waveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>(
			type, outputPrefix, param.timestep);

#ifdef USE_MPI
		m_waveFieldWriter->setComm(m_comm);
#endif // USE_MPI

		m_waveFieldWriter->init(variables, std::vector<const char*>(), true, true, true);
		m_waveFieldWriter->setMesh(
			info.bufferSize(param.bufferIds[CELLS]) / (4*sizeof(unsigned int)),
			static_cast<const unsigned int*>(info.buffer(param.bufferIds[CELLS])),
			info.bufferSize(param.bufferIds[VERTICES]) / (3*sizeof(double)),
			static_cast<const double*>(info.buffer(param.bufferIds[VERTICES])),
			param.timestep != 0);

		logInfo(rank) << "High order output initialized";

		//
		// Low order I/O
		//
		if (param.bufferIds[LOWCELLS] >= 0) {
			// Pstrain or Integrated quantities enabled
			m_lowOutputFlags = static_cast<const bool*>(info.buffer(param.bufferIds[LOW_OUTPUT_FLAGS]));
			// Variables
			std::vector<const char*> lowVariables;
			const char* lowVarNames[NUM_LOWVARIABLES] = {
				"ep_xx",
				"ep_yy",
				"ep_zz",
				"ep_xy",
				"ep_yz",
				"ep_xz",
				"eta",
				"int_sigma_xx",
				"int_sigma_yy",
				"int_sigma_zz",
				"int_sigma_xy",
				"int_sigma_yz",
				"int_sigma_xz",
				"displacement_x",
				"displacement_y",
				"displacement_z"
			};

			for (size_t i = 0; i < NUM_LOWVARIABLES; i++) {
				if (m_lowOutputFlags[i]) {
					lowVariables.push_back(lowVarNames[i]);
				}
			}

			m_lowWaveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON, double, real>(
				type, (std::string(outputPrefix)+"-low").c_str());

#ifdef USE_MPI
		m_lowWaveFieldWriter->setComm(m_comm);
#endif // USE_MPI

			m_lowWaveFieldWriter->init(lowVariables, std::vector<const char*>());
			m_lowWaveFieldWriter->setMesh(
				info.bufferSize(param.bufferIds[LOWCELLS]) / (4*sizeof(unsigned int)),
				static_cast<const unsigned int*>(info.buffer(param.bufferIds[LOWCELLS])),
				info.bufferSize(param.bufferIds[LOWVERTICES]) / (3*sizeof(double)),
				static_cast<const double*>(info.buffer(param.bufferIds[LOWVERTICES])),
				param.timestep != 0);

			logInfo(rank) << "Low order output initialized";
		}

		// Save ids for the variables
		m_variableBufferIds[0] = param.bufferIds[VARIABLE0];
		m_variableBufferIds[1] = param.bufferIds[LOWVARIABLE0];

		logInfo(rank) << "Initializing XDMF wave field output. Done.";
#ifdef USE_MPI
	}
	// End the if statement
#endif // USE_MPI
	}

	void setClusteringData(const unsigned *Clustering) {
	  m_waveFieldWriter->writeClusteringInfo(Clustering);
	}

	void exec(const async::ExecInfo &info, const WaveFieldParam &param)
	{
#ifdef USE_MPI
	// Execute this function only if m_waveFieldWriter is initialized
		if (m_waveFieldWriter != 0L) {
#endif // USE_MPI
		m_stopwatch.start();

		// High order output
		m_waveFieldWriter->addTimeStep(param.time);

		unsigned int nextId = 0;
		for (unsigned int i = 0; i < m_numVariables; i++) {
			if (m_outputFlags[i]) {
				m_waveFieldWriter->writeCellData(nextId,
					static_cast<const real*>(info.buffer(m_variableBufferIds[0]+nextId)));

				nextId++;
			}
		}

		m_waveFieldWriter->flush();

		// Low order output
		if (m_lowWaveFieldWriter) {
			m_lowWaveFieldWriter->addTimeStep(param.time);

		nextId = 0;
		for (unsigned int i = 0; i < NUM_LOWVARIABLES; i++) {
			if (m_lowOutputFlags[i]) {
				m_lowWaveFieldWriter->writeCellData(nextId,
					static_cast<const real*>(info.buffer(m_variableBufferIds[1]+nextId)));

			nextId++;
			}
		}

			m_lowWaveFieldWriter->flush();
		}

		m_stopwatch.pause();
#ifdef USE_MPI
		}
#endif // USE_MPI
	}

	void finalize()
	{
		if (m_waveFieldWriter) {
			m_stopwatch.printTime("Time wave field writer backend:"
#ifdef USE_MPI
				, m_comm
#endif // USE_MPI
			);
		}

#ifdef USE_MPI
		if (m_comm != MPI_COMM_NULL) {
			MPI_Comm_free(&m_comm);
			m_comm = MPI_COMM_NULL;
		}
#endif // USE_MPI

		delete m_waveFieldWriter;
		m_waveFieldWriter = 0L;
		delete m_lowWaveFieldWriter;
		m_lowWaveFieldWriter = 0L;
	}

public:
	static const unsigned int NUM_PLASTICITY_VARIABLES = 7;
	static const unsigned int NUM_INTEGRATED_VARIABLES = 9;
	static const unsigned int NUM_LOWVARIABLES = NUM_PLASTICITY_VARIABLES+NUM_INTEGRATED_VARIABLES;
};

}

}

#endif // WAVE_FIELD_WRITER_EXECUTOR_H

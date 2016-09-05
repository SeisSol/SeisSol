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

#ifndef WAVE_FIELD_WRITER_EXECUTOR_H
#define WAVE_FIELD_WRITER_EXECUTOR_H

#include "Parallel/MPI.h"

#include <cassert>
#include <vector>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "async/ExecInfo.h"

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
	LOWVARIABLE0,
	BUFFERTAG_MAX = LOWVARIABLE0
};

struct WaveFieldInitParam
{
	int timestep;

	int bufferIds[BUFFERTAG_MAX+1];
};

struct WaveFieldParam
{
	double time;
};

class WaveFieldWriterExecutor
{
private:
	/** The XMDF Writer used for the wave field */
	xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_waveFieldWriter;

	/** The XDMF Writer for low order data */
	xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_lowWaveFieldWriter;

	/** Buffer id for the first variable for high and low order output */
	unsigned int m_variableBufferIds[2];

	/** The total number of (high order) variables */
	unsigned int m_numVariables;

	/** Flag indicated which variables should be written */
	const bool* m_outputFlags;

#ifdef USE_MPI
	/** The MPI communicator for the XDMF writer */
	MPI_Comm m_comm;
#endif // USE_MPI

public:
	WaveFieldWriterExecutor()
		: m_waveFieldWriter(0L),
		  m_lowWaveFieldWriter(0L),
		  m_numVariables(0),
		  m_outputFlags(0L)
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

		const int rank = seissol::MPI::mpi.rank();

		const char* outputPrefix = static_cast<const char*>(info.buffer(param.bufferIds[OUTPUT_PREFIX]));

		//
		// High order I/O
		//
		m_numVariables = info.bufferSize(param.bufferIds[OUTPUT_FLAGS]) / sizeof(bool);
		m_outputFlags = static_cast<const bool*>(info.buffer(param.bufferIds[OUTPUT_FLAGS]));

		assert(m_numVariables <= 9);

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
			if (m_outputFlags[i])
				variables.push_back(varNames[i]);
		}

		// Initialize the I/O handler and write the mesh
		m_waveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
			rank, outputPrefix, variables, param.timestep);

#ifdef USE_MPI
		MPI_Comm_dup(seissol::MPI::mpi.comm(), &m_comm);
		m_waveFieldWriter->setComm(m_comm);
#endif // USE_MPI

		m_waveFieldWriter->init(
			info.bufferSize(param.bufferIds[CELLS]) / (4*sizeof(unsigned int)),
			static_cast<const unsigned int*>(info.buffer(param.bufferIds[CELLS])),
			info.bufferSize(param.bufferIds[VERTICES]) / (3*sizeof(double)),
			static_cast<const double*>(info.buffer(param.bufferIds[VERTICES])),
			true);

		logInfo(rank) << "High order output initialized";

		//
		// Low order I/O
		//
		if (param.bufferIds[LOWCELLS] >= 0) {
			// Pstrain enabled

			// Variables
			std::vector<const char*> lowVariables(NUM_LOWVARIABLES);
			lowVariables[0] = "ep_xx";
			lowVariables[1] = "ep_yy";
			lowVariables[2] = "ep_zz";
			lowVariables[3] = "ep_xy";
			lowVariables[4] = "ep_yz";
			lowVariables[5] = "ep_xz";
			lowVariables[6] = "eta";

			m_lowWaveFieldWriter = new xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>(
				rank, (std::string(outputPrefix)+"-low").c_str(), lowVariables, param.timestep);

#ifdef USE_MPI
			m_lowWaveFieldWriter->setComm(m_comm);
#endif // USE_MPI

			m_lowWaveFieldWriter->init(
				info.bufferSize(param.bufferIds[LOWCELLS]) / (4*sizeof(unsigned int)),
				static_cast<const unsigned int*>(info.buffer(param.bufferIds[LOWCELLS])),
				info.bufferSize(param.bufferIds[LOWVERTICES]) / (3*sizeof(double)),
				static_cast<const double*>(info.buffer(param.bufferIds[LOWVERTICES])),
				true);

			logInfo(rank) << "Low order output initialized";
		}

		// Save ids for the variables
		m_variableBufferIds[0] = param.bufferIds[VARIABLE0];
		m_variableBufferIds[1] = param.bufferIds[LOWVARIABLE0];

		logInfo(rank) << "Initializing XDMF wave field output. Done.";
	}

	void exec(const async::ExecInfo &info, const WaveFieldParam &param)
	{
		// High order output
		m_waveFieldWriter->addTimeStep(param.time);

		unsigned int nextId = 0;
		for (unsigned int i = 0; i < m_numVariables; i++) {
			if (m_outputFlags[i]) {
				m_waveFieldWriter->writeData(nextId,
					static_cast<const double*>(info.buffer(m_variableBufferIds[0]+nextId)));

				nextId++;
			}
		}

		m_waveFieldWriter->flush();

		// Low order output
		if (m_lowWaveFieldWriter) {
			m_lowWaveFieldWriter->addTimeStep(param.time);

			for (unsigned int i = 0; i < NUM_LOWVARIABLES; i++) {
				m_lowWaveFieldWriter->writeData(i,
					static_cast<const double*>(info.buffer(m_variableBufferIds[1]+i)));
			}

			m_lowWaveFieldWriter->flush();
		}
	}

	void finalize()
	{
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
	static const unsigned int NUM_LOWVARIABLES = 7;
};

}

}

#endif // WAVE_FIELD_WRITER_EXECUTOR_H
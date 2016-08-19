/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
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

#ifndef WAVE_FIELD_WRITER_H
#define WAVE_FIELD_WRITER_H

#include "Parallel/MPI.h"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include "utils/logger.h"

#include "xdmfwriter/XdmfWriter.h"

#include "async/Module.h"

#include "Geometry/refinement/VariableSubSampler.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol
{

namespace writer
{

/** Buffer tags for asynchronous IO */
enum BufferTags {
	OUTPUT_PREFIX,
	CELLS,
	VERTICES,
	VARIABLE0,
	VARIABLE1,
	VARIABLE2,
	VARIABLE3,
	VARIABLE4,
	VARIABLE5,
	VARIABLE6,
	VARIABLE7,
	VARIABLE8,
	LOWCELLS,
	LOWVERTICES,
	LOWVARIABLE0,
	LOWVARIABLE1,
	LOWVARIABLE2,
	LOWVARIABLE3,
	LOWVARIABLE4,
	LOWVARIABLE5,
	LOWVARIABLE6,
	BUFFERTAG_MAX = LOWVARIABLE6
};

const static unsigned int MAX_VARIABLES = VARIABLE8 - VARIABLE0 + 1;
const static unsigned int MAX_LOWVARIABLES = LOWVARIABLE6 - LOWVARIABLE0 + 1;

struct WaveFieldInitParam
{
	int timestep;
	unsigned int numVars;
#ifdef USE_ASYNC_MPI
	int bufferIds[BUFFERTAG_MAX+1];
#else // USE_ASYNC_MPI
	// Refined mesh structure
	size_t numCells;
	const unsigned int* cells;
	size_t numVertices;
	const double* vertices;

	// Original mesh structure
	size_t numLowCells;
	const unsigned int* lowCells;
	size_t numLowVertices;
	const double* lowVertices;
#endif // USE_ASYNC_MPI
};

struct WaveFieldParam
{
	double time;
};

class WaveFieldWriter : private async::Module<WaveFieldWriter, WaveFieldInitParam, WaveFieldParam>
{
	/** True if wave field output is enabled */
	bool m_enabled;

#ifdef USE_MPI
	/** The communicator for this output (a duplicate of MPI_COMM_WORLD) */
	MPI_Comm m_comm;
#endif // USE_MPI

	/** List of all buffer ids */
	int m_bufferIds[BUFFERTAG_MAX+1];

	/** The output prefix for the filename */
	std::string m_outputPrefix;

	/** The XMDF Writer used for the wave field */
	xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_waveFieldWriter;

	/** The XDMF Writer for low order data */
	xdmfwriter::XdmfWriter<xdmfwriter::TETRAHEDRON>* m_lowWaveFieldWriter;

	/** The variable subsampler for the refined mesh */
	refinement::VariableSubsampler<double>* m_variableSubsampler;

	/** Refined number of cells */
	unsigned int m_numCells;

	/** Unrefined (low order) number of cells */
	unsigned int m_numLowCells;

	/** Number of variables */
	unsigned int m_numVariables;

	/** Pointer to the degrees of freedom */
	const double* m_dofs;

	/** Pointer to the plastic strain */
	const double* m_pstrain;

	/** Mapping from the cell order to dofs order */
	const unsigned int* m_map;

	/** Time of the last output (makes sure output is not written twice at the end) */
	double m_lastTimeStep;

	/** The tolerance in the time for ignoring duplicate time steps */
	double m_timeTolerance;

	/** The current output time step */
	unsigned int m_timestep;

	/** Buffer required to extract the output data from the unknowns */
	double* m_outputBuffer;

	/** Flag indicated which variables should be written */
	std::vector<bool> m_outputFlags;

public:
	WaveFieldWriter()
		: m_enabled(false),
#ifdef USE_MPI
		  m_comm(MPI_COMM_NULL),
#endif // USE_MPI
		  m_waveFieldWriter(0L), m_lowWaveFieldWriter(0L),
		  m_variableSubsampler(0L),
		  m_numCells(0), m_numLowCells(0),
		  m_numVariables(0),
		  m_dofs(0L), m_pstrain(0L),
		  m_map(0L),
		  m_lastTimeStep(-1),
		  m_timeTolerance(0),
		  m_timestep(0),
		  m_outputBuffer(0L)
	{
		std::fill(m_bufferIds, m_bufferIds+BUFFERTAG_MAX+1, -1);
	}

	/**
	 * Activate the wave field output
	 */
	void enable()
	{
		m_enabled = true;
	}

	/**
	 * @return True if wave field output is enabled, false otherwise
	 */
	bool isEnabled() const
	{
		return m_enabled;
	}

	/**
	 * Set the output prefix for the filename
	 */
	void setFilename(const char* outputPrefix)
	{
		m_outputPrefix = outputPrefix;
	}

	/**
	 * This is called on the executor instead of {@link init()}
	 */
	void setUp()
	{
		setExecutor(*this);
	}

	/**
	 * Initialize the wave field ouput
	 *
	 * @param map The mapping from the cell order to dofs order
	 * @param timeTolerance The tolerance in the time for ignoring duplicate time steps
	 */
	void init(unsigned int numVars, int order, int numAlignedDOF,
			const MeshReader &meshReader,
			const double* dofs,  const double* pstrain,
			const unsigned int* map,
			int refinement, int timestep,
			double timeTolerance);

	/**
	 * Initialization on the executor
	 */
	void execInit(const WaveFieldInitParam &param);

	/**
	 * @return The current time step of the wave field output
	 */
	unsigned int timestep() const
	{
		return m_timestep;
	}

	/**
	 * Write a time step
	 */
	void write(double time)
	{
		SCOREP_USER_REGION("WaveFieldWriter_write", SCOREP_USER_REGION_TYPE_FUNCTION);

		if (!m_enabled)
			return;

		const int rank = seissol::MPI::mpi.rank();

		if (time <= m_lastTimeStep + m_timeTolerance) {
			// Ignore duplicate time steps. Might happen at the end of a simulation
			logInfo(rank) << "Ignoring duplicate time step at time " << time;
			return;
		}

		SCOREP_USER_REGION_DEFINE(r_wait);
		SCOREP_USER_REGION_BEGIN(r_wait, "wavfieldwriter_wait", SCOREP_USER_REGION_TYPE_COMMON);
		logInfo(rank) << "Waiting for last wave field.";
		wait();
		SCOREP_USER_REGION_END(r_wait);

		logInfo(rank) << "Writing wave field at time" << utils::nospace <<  time << '.';

#if !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))
		// We can use the exec function in synchronous mode since we have only
		// one buffer avaiable
        m_waveFieldWriter->addTimeStep(time);
#endif // !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))

		for (unsigned int i = 0; i < m_numVariables; i++) {
			if (!m_outputFlags[i])
				continue;

			m_variableSubsampler->get(m_dofs, m_map, i, m_outputBuffer);

#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
			assert(m_bufferIds[VARIABLE0+i] >= 0);
			fillBuffer(m_bufferIds[VARIABLE0+i], m_outputBuffer, m_numCells*sizeof(double));
#else // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
			m_waveFieldWriter->writeData(i, m_outputBuffer);
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
		}

#if !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))
		m_waveFieldWriter->flush();
#endif // !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))

		if (m_pstrain) {
#if !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))
			m_lowWaveFieldWriter->addTimeStep(time);
#endif // !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))

			for (unsigned int i = 0; i < MAX_LOWVARIABLES; i++) {
#ifdef _OPENMP
				#pragma omp parallel for schedule(static)
#endif // _OPENMP
				for (unsigned int j = 0; j < m_numLowCells; j++)
					m_outputBuffer[j] = m_pstrain[m_map[j] * 7 + i];

#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
				assert(m_bufferIds[LOWVARIABLE0+i] >= 0);
				fillBuffer(m_bufferIds[LOWVARIABLE0+i], m_outputBuffer, m_numLowCells*sizeof(double));
#else // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
				m_lowWaveFieldWriter->writeData(i, m_outputBuffer);
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
			}

#if !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))
			m_lowWaveFieldWriter->flush();
#endif // !(defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD))
		}

#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
		WaveFieldParam param;
		param.time = time;
		call(param);
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)

		// Update last time step
		m_lastTimeStep = time;
		m_timestep++;

		logInfo(rank) << "Writing wave field at time" << utils::nospace << time << ". Done.";
	}

	void exec(const WaveFieldParam& param)
	{
		// High order output
		m_waveFieldWriter->addTimeStep(param.time);

		for (unsigned int i = 0; i < m_numVariables; i++) {
			if (m_bufferIds[VARIABLE0+i] >= 0) {
				assert(bufferSize(m_bufferIds[VARIABLE0+i]) > 0);
				m_waveFieldWriter->writeData(i,
						static_cast<const double*>(buffer(m_bufferIds[VARIABLE0+i])));
			}
		}

		m_waveFieldWriter->flush();

		// Low order output
		if (m_lowWaveFieldWriter) {
			m_lowWaveFieldWriter->addTimeStep(param.time);

			for (unsigned int i = 0; i < 7; i++) {
				if (m_bufferIds[LOWVARIABLE0+i] >= 0) {
					assert(bufferSize(m_bufferIds[LOWVARIABLE0+i]) > 0);
					m_lowWaveFieldWriter->writeData(i,
							static_cast<const double*>(buffer(m_bufferIds[LOWVARIABLE0+i])));
				}
			}

			m_lowWaveFieldWriter->flush();
		}
	}

	/**
	 * Close wave field writer and free resources
	 */
	void close()
	{
		// Cleanup the executor
		wait();
		finalize();

		if (!m_enabled)
			return;

		delete m_variableSubsampler;
		m_variableSubsampler = 0L;
		delete m_outputBuffer;
		m_outputBuffer = 0L;
	}

	void tearDown()
	{
		if (!m_enabled)
			return;

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
};

}

}

#endif // WAVE_FIELD_WRITER_H

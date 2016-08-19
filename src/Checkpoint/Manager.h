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
 * Manages the checkpoints
 */

#ifndef CHECKPOINT_MANAGER_H
#define CHECKPOINT_MANAGER_H

#include "Parallel/MPI.h"

#include <cstring>
#include <string>

#include "async/Module.h"

#include "Wavefield.h"
#include "Fault.h"
#include "posix/Fault.h"
#include "posix/Wavefield.h"
#include "h5/Wavefield.h"
#include "h5/Fault.h"
#include "mpio/Wavefield.h"
#include "mpio/WavefieldAsync.h"
#include "mpio/Fault.h"
#include "mpio/FaultAsync.h"
#ifdef USE_SIONLIB
#include "sionlib/Fault.h"
#include "sionlib/Wavefield.h"
#endif // USE_SIONLIB

namespace seissol
{

namespace checkpoint
{

/** Checkpoint backend types */
enum Backend {
	POSIX,
	HDF5,
	MPIO,
	MPIO_ASYNC,
	SIONLIB,
	DISABLED
};

/** Buffer tags for asynchronous IO */
enum BufferTags {
	FILENAME,
	DOFS,
	MU,
	SLIP_RATE1,
	SLIP_RATE2,
	SLIP,
	SLIP1,
	SLIP2,
	STATE,
	STRENGTH,
	BUFFERTAG_MAX = STRENGTH
};

/**
 * Initialization parameters for checkpoints
 */
struct CheckpointInitParam
{
#ifdef USE_ASYNC_MPI
	int bufferIds[BUFFERTAG_MAX+1];
	Backend backend;
#endif // USE_ASYNC_MPI
	unsigned int numBndGP;
	bool loaded;
};

/**
 * Parameters for checkpoints
 */
struct CheckpointParam
{
	double time;
	int waveFieldTimeStep;
	int faultTimeStep;
};

class Manager : private async::Module<Manager, CheckpointInitParam, CheckpointParam>
{
private:
	/** The backend that should be used */
	Backend m_backend;

	/** The filename for the checkpoints */
	std::string m_filename;

	/** The wave field checkpoint */
	Wavefield *m_waveField;

	/** The dynamic rupture checkpoint */
	Fault *m_fault;

	/** List of all buffer ids */
	int m_bufferIds[BUFFERTAG_MAX+1];

	/** Number of DOFs */
	unsigned int m_numDofs;

	/** Number of DR DOFs */
	unsigned int m_numDRDofs;

	/** Pointer to the dofs (required for async IO) */
	const real* m_dofs;

	/** Pointer to the DR dofs (required for async IO */
	const double* m_drDofs[8];

public:
	Manager()
		: m_backend(DISABLED),
		  m_waveField(0L), m_fault(0L),
		  m_numDofs(0), m_numDRDofs(0),
		  m_dofs(0L)
	{
		std::fill(m_bufferIds, m_bufferIds+BUFFERTAG_MAX, -1);
		std::fill(m_drDofs, m_drDofs+8, static_cast<const double*>(0L));
	}

	virtual ~Manager()
	{
		delete m_waveField;
		delete m_fault;
	}

	void setBackend(Backend backend)
	{
		m_backend = backend;
	}

	/**
	 * Set the filename prefix for checkpointing
	 *
	 * @param filename The filename prefix
	 */
	void setFilename(const char* filename)
	{
		m_filename = filename;
	}


	/**
	 * This is called on the executor instead of {@link init()}
	 */
	void setUp()
	{
		setExecutor(*this);
	}

	/**
	 * Initialize checkpointing and load the last checkpoint if present
	 *
	 * @return True is a checkpoint was loaded, false otherwise
	 */
	bool init(real* dofs, unsigned int numDofs,
			double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1, double* slip2,
			double* state, double* strength, unsigned int numSides, unsigned int numBndGP,
			double &time, int &waveFieldTimeStep, int &faultTimeStep);

	/**
	 * Initializer on the executor
	 */
	void execInit(const CheckpointInitParam &param)
	{
#ifdef USE_ASYNC_MPI
		memcpy(m_bufferIds, param.bufferIds, sizeof(m_bufferIds));
		const char* filename = static_cast<const char*>(buffer(m_bufferIds[FILENAME]));
		m_backend = param.backend;

		initBackend();

		m_waveField->setFilename(filename);
		m_fault->setFilename(filename);

		m_waveField->init(bufferSize(m_bufferIds[DOFS]) / sizeof(real));
		m_fault->init(bufferSize(m_bufferIds[MU]) / param.numBndGP / sizeof(double), param.numBndGP);
#endif // USE_ASYNC_MPI

		if (param.loaded) {
			m_waveField->setLoaded();
			m_fault->setLoaded();
		}

		const real* dofs;
		const double* drDofs[8];
#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
		dofs = static_cast<const real*>(buffer(m_bufferIds[DOFS]));
		for (unsigned int i = 0; i < 8; i++)
			drDofs[i] = static_cast<const double*>(buffer(m_bufferIds[MU + i]));
#else // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
		dofs = m_dofs;
		memcpy(drDofs, m_drDofs, 8*sizeof(double*));
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)

		m_waveField->initLate(dofs);
		m_fault->initLate(drDofs[0], drDofs[1], drDofs[2], drDofs[3], drDofs[4], drDofs[5],
			drDofs[6], drDofs[7]);
	}

	/**
	 * Write a checkpoint for the current time
	 *
	 * @param time The current time
	 * @param waveFieldTimeStep The time step of the wave field writer
	 * @param faultTimeStep The time step of the fault writer
	 */
	void write(double time, int waveFieldTimeStep, int faultTimeStep)
	{
		SCOREP_USER_REGION("CheckpointManager_write", SCOREP_USER_REGION_TYPE_FUNCTION);

		if (m_backend == DISABLED)
			return;

		const int rank = seissol::MPI::mpi.rank();

		SCOREP_USER_REGION_DEFINE(r_wait);
		SCOREP_USER_REGION_BEGIN(r_wait, "checkpointmanager_wait", SCOREP_USER_REGION_TYPE_COMMON);
		logInfo(rank) << "Waiting for last checkpoint.";
		wait();
		SCOREP_USER_REGION_END(r_wait);

		logInfo(rank) << "Writing checkpoint at time" << utils::nospace <<  time << '.';

#if defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)
		// Fill buffers
		fillBuffer(m_bufferIds[DOFS], m_dofs, m_numDofs * sizeof(real));
		for (unsigned int i = 0; i < 8; i++)
			fillBuffer(m_bufferIds[MU+i], m_drDofs[i], m_numDRDofs * sizeof(double));
#endif // defined(USE_ASYNC_MPI) || defined(USE_ASYNC_THREAD)

		CheckpointParam param;
		param.time = time;
		param.waveFieldTimeStep = waveFieldTimeStep;
		param.faultTimeStep = faultTimeStep;
		call(param);

		logInfo(rank) << "Writing checkpoint at time" << utils::nospace << time << ". Done.";
	}

	/**
	 * Called on the executor
	 */
	void exec(const CheckpointParam &param)
	{
		m_waveField->write(param.time, param.waveFieldTimeStep);
		m_fault->write(param.faultTimeStep);

		// Update both links at the "same" time
		m_waveField->updateLink();
		m_fault->updateLink();

		// Prepare next checkpoint (only for async checkpoints)
		m_waveField->writePrepare(param.time, param.waveFieldTimeStep);
		m_fault->writePrepare(param.faultTimeStep);
	}

	/**
	 * Close checkpointing
	 */
	void close()
	{
		// Cleanup/terminate the executor
		wait();

		if (m_backend != DISABLED) {
			m_waveField->close();
			m_fault->close();
		}

		finalize(); // Will also call tearDown for the local process
	}

	void tearDown()
	{
		if (m_backend == DISABLED)
			return;

		delete m_waveField;
		m_waveField = 0L;
		delete m_fault;
		m_fault = 0L;
	}

private:
	/**
	 * Initialize the backend
	 */
	void initBackend()
	{
		switch (m_backend) {
		case POSIX:
			m_waveField = new posix::Wavefield();
			m_fault = new posix::Fault();
			break;
		case HDF5:
			m_waveField = new h5::Wavefield();
			m_fault = new h5::Fault();
			break;
		case MPIO:
			m_waveField = new mpio::Wavefield();
			m_fault = new mpio::Fault();
			break;
		case MPIO_ASYNC:
			m_waveField = new mpio::WavefieldAsync();
			m_fault = new mpio::FaultAsync();
			break;
		case SIONLIB:
#ifdef USE_SIONLIB
			m_waveField = new sionlib::Wavefield();
			m_fault = new sionlib::Fault();
			break;
#else //USE_SIONLIB
			logError() << "SIONlib checkpoint backend unsupported";
			break;
#endif //USE_SIONLIB
		default:
			logError() << "Unsupported checkpoint backend";
		}
	}
};

}

}

#endif // CHECKPOINT_MANAGER_H

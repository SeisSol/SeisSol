/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

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
#endif
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

class Manager
{
private:
	/** The wave field checkpoint */
	Wavefield *m_waveField;

	/** The dynamic rupture checkpoint */
	Fault *m_fault;

public:
	Manager()
		: m_waveField(0L), m_fault(0L)
	{
	}

	virtual ~Manager()
	{
		delete m_waveField;
		delete m_fault;
	}

	/**
	 * Enables checkpointing
	 */
	void setBackend(Backend backend)
	{
		switch (backend) {
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
			m_waveField = new posix::Wavefield();
			m_fault = new posix::Fault();
			logError() << "sionlib checkpoint backend unsupported - using posix instead";
			break;
#endif //USE_SIONLIB
		default:
			logError() << "Unsupported checkpoint backend";
		}
	}

	/**
	 * Set the filename prefix for checkpointing
	 *
	 * @param filename The filename prefix
	 */
	void setFilename(const char* filename)
	{
		if (!m_waveField)
			return;

		m_waveField->setFilename(filename);
		m_fault->setFilename(filename);
	}

	/**
	 * Initialize checkpointing and load the last checkpoint if present
	 *
	 * @return True is a checkpoint was loaded, false otherwise
	 */
	bool init(real* dofs, unsigned int numDofs,
			double* mu, double* slipRate1, double* slipRate2, double* slip, double* slip1, double* slip2,

			double* state, double* strength, unsigned int numSides, unsigned int numBndGP,
			double &time, int &waveFieldTimeStep, int &faultTimeStep)
	{
		if (!m_waveField)
			return false;

		int exists = m_waveField->init(dofs, numDofs);
		exists &= m_fault->init(mu, slipRate1, slipRate2, slip, slip1, slip2,
				state, strength, numSides, numBndGP);

		// Make sure all rank think the same about the existing checkpoint
#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &exists, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
#endif // USE_MPI

		// Load checkpoint?
		if (exists) {
			m_waveField->load(time, waveFieldTimeStep);
			m_fault->load(faultTimeStep);
		}

		m_waveField->initLate();
		m_fault->initLate();

		return exists;
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
		if (!m_waveField)
			return;

		m_waveField->write(time, waveFieldTimeStep);
		m_fault->write(faultTimeStep);

		// Update both links at the "same" time
		m_waveField->updateLink();
		m_fault->updateLink();

		// Prepare next checkpoint (only for async checkpoints)
		m_waveField->writePrepare(time, waveFieldTimeStep);
		m_fault->writePrepare(faultTimeStep);
	}

	/**
	 * Close checkpointing
	 */
	void close()
	{
		if (!m_waveField)
			return;

		m_waveField->close();
		m_fault->close();
	}
};

}

}

#endif // CHECKPOINT_MANAGER_H

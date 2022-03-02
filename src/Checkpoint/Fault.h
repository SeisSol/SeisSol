/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2017, SeisSol Group
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
 * Common interface for fault checkpoints
 */

#ifndef CHECKPOINT_FAULT_H
#define CHECKPOINT_FAULT_H

#include "Parallel/MPI.h"

#include <cassert>
#include <string>

#include "CheckPoint.h"

namespace seissol
{

namespace checkpoint
{

/**
 * Common interface for fault checkpoints
 */
class Fault : virtual public CheckPoint
{
protected:
	static const unsigned int NUM_VARIABLES = 8;

private:
	/** Pointers to fault data */
	const double* m_data[NUM_VARIABLES];

	/** Number of dynamic rupture sides on this rank */
	unsigned int m_numSides;

	/** Number of boundary points per side */
	unsigned int m_numBndGP;

public:
	Fault(unsigned long identifier)
		: CheckPoint(identifier),
		  m_numSides(0), m_numBndGP(0)
	{}

	virtual ~Fault() {}

	/**
	 * @return True of a valid checkpoint is available
	 */
	virtual bool init(unsigned int numSides, unsigned int numBndGP,
			unsigned int groupSize = 1)
	{
		const int rank = seissol::MPI::mpi.rank();

		logInfo(rank) << "Initializing fault check pointing";

#ifdef USE_MPI
		// Get the communicator (must be called by all ranks)
		MPI_Comm comm;
		MPI_Comm_split(seissol::MPI::mpi.comm(), (numSides == 0 ? MPI_UNDEFINED : 0), rank, &comm);
#endif // USE_MPI

		if (numSides == 0)
			return true;

#ifdef USE_MPI
		setComm(comm);
#endif // USE_MPI

		// Compute sum and offset for this group
		setGroupSumOffset(numSides, groupSize);

		m_numSides = numSides;
		m_numBndGP = numBndGP;

		return false;
	}

	/**
	 * @param[out] timestepFault Time step of the fault writer in the checkpoint
	 *  (if the fault writer was active)
	 */
	virtual void load(int &timestepFault, double* mu, double* slipRate1, double* slipRate2,
		double* slip, double* slip1, double* slip2, double* state, double* strength) = 0;

	/**
	 * @copydoc CheckPoint::initLate
	 */
	virtual void initLate(const double* mu, const double* slipRate1, const double* slipRate2,
		const double* slip, const double* slip1, const double* slip2, const double* state, const double* strength)
	{
		if (numSides() == 0)
			return;

		m_data[0] = mu;
		m_data[1] = slipRate1;
		m_data[2] = slipRate2;
		m_data[3] = slip;
		m_data[4] = slip1;
		m_data[5] = slip2;
		m_data[6] = state;
		m_data[7] = strength;

		createFiles();
	}

	/**
	 * Prepare writing a checkpoint
	 *
	 * Overwrite this function to start an asynchronous checkpoint write
	 */
	virtual void writePrepare(int timestepFault)
	{
	}

	virtual void write(int timestepFault) = 0;

	void updateLink()
	{
		if (m_numSides == 0)
			return;

		CheckPoint::updateLink();
	}

protected:
	void initFilename(const char* filename, const char* extension)
	{
		std::string file(filename);
		std::string ext;
		if (extension)
			ext = std::string(".") + extension;
		if (utils::StringUtils::endsWith(file, ext))
			utils::StringUtils::replaceLast(file, ext, "-fault");
		else
			file += "-fault";

		CheckPoint::initFilename(file.c_str(), extension);
	}

	virtual const char* fname() const
	{
		return "cp-f";
	}

	void setData(const double* mu, const double* slipRate1, const double* slipRate2,
			const double* slip, const double* slip1, const double* slip2,
			const double* state, const double* strength)
	{
		m_data[0] = mu;
		m_data[1] = slipRate1;
		m_data[2] = slipRate2;
		m_data[3] = slip;
		m_data[4] = slip1;
		m_data[5] = slip2;
		m_data[6] = state;
		m_data[7] = strength;
	}

	const double* data(unsigned int var) const
	{
		assert(var < NUM_VARIABLES);
		return m_data[var];
	}

	unsigned int numSides() const
	{
		return m_numSides;
	}

	unsigned int numBndGP() const
	{
		return m_numBndGP;
	}

	/** Names of the different variables we need to store */
	static const char* VAR_NAMES[NUM_VARIABLES];
};

}

}

#endif // CHECKPOINT_FAULT_H

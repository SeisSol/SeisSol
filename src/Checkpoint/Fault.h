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
 * Common interface for fault checkpoints
 */

#ifndef CHECKPOINT_FAULT_H
#define CHECKPOINT_FAULT_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

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
	static const int NUM_VARIABLES = 6;

private:
	/** Pointers to fault data */
	double* m_data[NUM_VARIABLES];

	/** Number of dynamic rupture sides on this rank */
	unsigned int m_numSides;

	/** Number of boundary points per side */
	unsigned int m_numBndGP;

public:
	Fault()
		: m_numSides(), m_numBndGP(0)
	{}

	virtual ~Fault() {}

	/**
	 * @return True of a valid checkpoint is available
	 */
	virtual bool init(double* mu, double* slipRate1, double* slipRate2, double* slip,
			double* state, double* strength,
			unsigned int numSides, unsigned int numBndGP)
	{
		int rank = 0;
#ifdef USE_MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif // USE_MPI

		logInfo(rank) << "Initializing fault check pointing";

#ifdef USE_MPI
		// Get the communicator (must be called by all ranks)
		MPI_Comm comm;
		MPI_Comm_split(MPI_COMM_WORLD, (numSides == 0 ? MPI_UNDEFINED : 0), rank, &comm);
#endif // USE_MPI

		if (numSides == 0)
			return true;

#ifdef USE_MPI
		setComm(comm);
#endif // USE_MPI

		// Save data pointers
		m_data[0] = mu;
		m_data[1] = slipRate1;
		m_data[2] = slipRate2;
		m_data[3] = slip;
		m_data[4] = state;
		m_data[5] = strength;
		m_numSides = numSides;
		m_numBndGP = numBndGP;

		return false;
	}

	/**
	 * @param[out] timestepFault Time step of the fault writer in the checkpoint
	 *  (if the fault writer was active)
	 */
	virtual void load(int &timestepFault) = 0;

	virtual void write(int timestepFault) = 0;

protected:
	void initFilename(const char* filename, const char* extension)
	{
		std::string file(filename);
		std::string ext = "." + std::string(extension);
		if (utils::StringUtils::endsWith(file, ext))
			utils::StringUtils::replaceLast(file, ext, "-fault");
		else
			file += "-fault";

		CheckPoint::initFilename(file.c_str(), extension);
	}

	double* data(unsigned int var)
	{
		assert(var < NUM_VARIABLES);
		return m_data[var];
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

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
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

#ifndef CHECKPOINT_H5_FAULT_H
#define CHECKPOINT_H5_FAULT_H

#ifndef USE_HDF
#include "Checkpoint/FaultDummy.h"
#else // USE_HDF

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <cstdio>
#include <string>

#include <hdf5.h>

#include "utils/arrayutils.h"
#include "utils/logger.h"

#include "CheckPoint.h"
#include "Checkpoint/Fault.h"
#include "Initializer/preProcessorMacros.fpp"
#include "Initializer/typedefs.hpp"

#endif // USE_HDF

namespace seissol
{

namespace checkpoint
{

namespace h5
{

#ifndef USE_HDF
typedef FaultDummy Fault;
#else // USE_HDF

class Fault : public CheckPoint, virtual public seissol::checkpoint::Fault
{
private:
	/** Pointer to fault data */
	double* m_mu;

	/** Pointer to fault data */
	double* m_slipRate1;

	/** Pointer to fault data */
	double* m_slipRate2;

	/** Pointer to fault data */
	double* m_slip;

	/** Pointer to fault data */
	double* m_state;

	/** Pointer to fault data */
	double* m_strength;

	/** Number of dynamic rupture sides on this rank */
	unsigned int m_numSides;

	/** Number of boundary points per side */
	unsigned int m_numBndGP;

	/** Identifiers of the HDF5 fault attributes */
	hid_t m_h5timestepFault[2];

	/**
	 * Identifiers of the main data set in the files
	 * @todo Make the size depend on the size of VAR_NAMES
	 */
	hid_t m_h5data[2][6];

	/** Identifiers for the file space of the data set */
	hid_t m_h5fSpaceData;

public:
	Fault()
		: m_mu(0L), m_slipRate1(0L), m_slipRate2(0L),
		  m_slip(0L), m_state(0L), m_strength(0L),
		  m_numSides(0), m_numBndGP(0),
		  m_h5fSpaceData(-1)
	{}

	void setFilename(const char* filename);

	bool init(double* mu, double* slipRate1, double* slipRate2, double* slip,
			double* state, double* strength,
			unsigned int numSides, unsigned int numBndGP);

	void initLate()
	{
		if (m_numSides == 0)
			return;

		CheckPoint::initLate();
	}

	/**
	 * @param[out] timestepFault Time step of the fault writer in the checkpoint
	 *  (if the fault writer was active)
	 */
	void load(int &timestepFault);

	void write(int timestepFault);

	void close()
	{
		if (m_numSides == 0)
			return;

		for (unsigned int i = 0; i < 2; i++) {
			checkH5Err(H5Aclose(m_h5timestepFault[i]));
			for (unsigned int j = 0; j < utils::ArrayUtils::size(VAR_NAMES); j++)
				checkH5Err(H5Dclose(m_h5data[i][j]));
		}
		checkH5Err(H5Sclose(m_h5fSpaceData));

		CheckPoint::close();
	}


private:
	bool validate(hid_t h5file) const;

	hid_t create(int odd, const char* filename);

	/** Names of the different variables we need to store */
	static const char* VAR_NAMES[6];
};

#endif // USE_HDF

}

}

}

#endif // CHECKPOINT_H5_FAULT_H

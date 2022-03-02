/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2017, SeisSol Group
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
	/** Identifiers of the HDF5 fault attributes */
	hid_t m_h5timestepFault[2];

	/**
	 * Identifiers of the main data set in the files
	 * @todo Make the size depend on the size of VAR_NAMES
	 */
	hid_t m_h5data[2][NUM_VARIABLES];

	/** Identifiers for the file space of the data set */
	hid_t m_h5fSpaceData;

public:
	Fault()
		: seissol::checkpoint::CheckPoint(IDENTIFIER),
		seissol::checkpoint::Fault(IDENTIFIER),
		CheckPoint(IDENTIFIER),
		m_h5fSpaceData(-1)
	{
		m_h5timestepFault[0] = m_h5timestepFault[1] = -1;
		for (unsigned int i = 0; i < NUM_VARIABLES; i++)
			m_h5data[0][i] = m_h5data[1][i] = -1;
	}

	bool init(unsigned int numSides, unsigned int numBndGP,
		unsigned int groupSize = 1);

	/**
	 * @param[out] timestepFault Time step of the fault writer in the checkpoint
	 *  (if the fault writer was active)
	 */
	void load(int &timestepFault, double* mu, double* slipRate1, double* slipRate2,
		double* slip, double* slip1, double* slip2, double* state, double* strength);

	void write(int timestepFault);

	void close()
	{
		if (numSides() == 0)
			return;

		if (m_h5timestepFault[0] >= 0) {
			for (unsigned int i = 0; i < 2; i++) {
				checkH5Err(H5Aclose(m_h5timestepFault[i]));
				for (unsigned int j = 0; j < NUM_VARIABLES; j++)
					checkH5Err(H5Dclose(m_h5data[i][j]));
			}
		}
		if (m_h5fSpaceData >= 0)
			checkH5Err(H5Sclose(m_h5fSpaceData));

		CheckPoint::close();
	}


private:
	bool validate(hid_t h5file) const;

	hid_t initFile(int odd, const char* filename);

private:
	static const unsigned long IDENTIFIER = 0x7A127;
};

#endif // USE_HDF

}

}

}

#endif // CHECKPOINT_H5_FAULT_H

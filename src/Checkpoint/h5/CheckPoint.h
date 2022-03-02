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
 * Main class for checkpoints
 */

#ifndef CHECKPOINT_H5_CHECK_POINT_H
#define CHECKPOINT_H5_CHECK_POINT_H

#include <string>

#include <hdf5.h>

#include "utils/logger.h"
#include "utils/path.h"

#include "H5ErrHandler.h"
#include "Checkpoint/CheckPoint.h"
#ifdef USE_MPI
#include "Checkpoint/MPIInfo.h"
#endif // USE_MPI
#include "Initializer/preProcessorMacros.fpp"

namespace seissol
{

namespace checkpoint
{

namespace h5
{

class CheckPoint : virtual public seissol::checkpoint::CheckPoint
{
private:
	/** Identifiers of the HDF5 files */
	hid_t m_h5files[2];

	/** Property list for data access */
	hid_t m_h5XferList;

public:
	CheckPoint(unsigned long identifier)
		: seissol::checkpoint::CheckPoint(identifier),
		m_h5XferList(-1)
	{
		m_h5files[0] = m_h5files[1] = -1;
	}

	virtual ~CheckPoint()
	{ }

	void setFilename(const char* filename)
	{
		initFilename(filename, "h5");
	}

	void close()
	{
		if (m_h5files[0] >= 0) {
			for (unsigned int i = 0; i < 2; i++)
				checkH5Err(H5Fclose(m_h5files[i]));
		}

		if (m_h5XferList >= 0)
			checkH5Err(H5Pclose(m_h5XferList));
	}

protected:
	/**
	 * Sets up transfer list for collective read/writes
	 */
	void setupXferList()
	{
		// Initialize HDF5 property lists
	#ifdef USE_MPI
		m_h5XferList = H5Pcreate(H5P_DATASET_XFER);
		checkH5Err(m_h5XferList);
		H5Pset_dxpl_mpio(m_h5XferList, H5FD_MPIO_COLLECTIVE);
	#else // USE_MPI
		m_h5XferList = H5P_DEFAULT;
	#endif // USE_MPI
	}

	bool exists()
	{
		if (!seissol::checkpoint::CheckPoint::exists())
			return false;

		hid_t h5file = open(linkFile());
		if (h5file < 0) {
			logWarning(rank()) << "Failed to open checkpoint file.";
			return false;
		}

		bool hasCheckpoint = validate(h5file);
		checkH5Err(H5Fclose(h5file));

		return hasCheckpoint;
	}

	void createFiles()
	{
		// Backup checkpoints from last run
		seissol::checkpoint::CheckPoint::createFiles();

		for (unsigned int i = 0; i < 2; i++) {
			m_h5files[i] = initFile(i, dataFile(i).c_str());

			// Sync file (required for performance measure)
			checkH5Err(H5Fflush(m_h5files[i], H5F_SCOPE_GLOBAL));
		}
	}

	/**
	 * Finalize checkpoint writing: Flush the file
	 */
	void finalizeCheckpoint()
	{
		EPIK_USER_REG(r_flush, "checkpoint_flush");
		SCOREP_USER_REGION_DEFINE(r_flush);
		EPIK_USER_START(r_flush);
		SCOREP_USER_REGION_BEGIN(r_flush, "checkpoint_flush", SCOREP_USER_REGION_TYPE_COMMON);

		checkH5Err(H5Fflush(m_h5files[odd()], H5F_SCOPE_GLOBAL));

		EPIK_USER_END(r_flush);
		SCOREP_USER_REGION_END(r_flush);
	}

	hid_t h5XferList() const
	{
		return m_h5XferList;
	}

	/**
	 * Open a check point file
	 *
	 * @return H5 identifier or -1 on failure
	 */
	hid_t open(const char* file, bool readonly = true)
	{
		hid_t h5plist = H5P_DEFAULT;
#ifdef USE_MPI
		h5plist = H5Pcreate(H5P_FILE_ACCESS);
		checkH5Err(h5plist);

		MPIInfo info;
		checkH5Err(H5Pset_fapl_mpio(h5plist, comm(), info.get()));
#endif // USE_MPI

		// Turn of error printing
		H5ErrHandler errHandler;

		// Do not check the id since this may file if no checkpoint exists
		hid_t h5file = H5Fopen(file, (readonly ? H5F_ACC_RDONLY : H5F_ACC_RDWR), h5plist);

		// Restore previous error handler
		errHandler.outputOn();

		checkH5Err(H5Pclose(h5plist));

		return h5file;
	}

	/**
	 * Validate an existing check point file
	 */
	virtual bool validate(hid_t h5file) const = 0;

	/**
	 * Create or open a check a checkpoint file for writing
	 *
	 * @param odd 1 if odd file should be created, 0 otherwise
	 * @param filename Path the the HDF5 file that should be created
	 * @return The id of the created file
	 */
	virtual hid_t initFile(int odd, const char* filename) = 0;

protected:
	template<typename T>
	static void checkH5Err(T status)
	{
		if (status < 0)
			logError() << "An error in the HDF5 checkpoint module occurred";
	}
};

}

}

}

#endif // CHECKPOINT_H5_CHECK_POINT_H

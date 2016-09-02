/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2014-2016, SeisSol Group
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

#include "Fault.h"

#ifdef USE_MPI
#include "Checkpoint/MPIInfo.h"
#endif // USE_MPI

bool seissol::checkpoint::h5::Fault::init(unsigned int numSides, unsigned int numBndGP,
		unsigned int groupSize)
{
	seissol::checkpoint::Fault::init(numSides, numBndGP);

	if (numSides == 0)
		return true;

	// Compute total number of cells and local offset
	setSumOffset(numSides);

	// Dataspace for the file
	hsize_t fileSize[2] = {numTotalElems(), numBndGP};
	m_h5fSpaceData = H5Screate_simple(2, fileSize, 0L);
	checkH5Err(m_h5fSpaceData);

	setupXferList();

	return exists();
}

void seissol::checkpoint::h5::Fault::load(int &timestepFault, double* mu, double* slipRate1, double* slipRate2,
	double* slip, double* slip1, double* slip2, double* state, double* strength)
{
	if (numSides() == 0)
		return;

	logInfo(rank()) << "Loading fault checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	hid_t h5file = open(linkFile());
	checkH5Err(h5file);

	// Attributes
	hid_t h5attr = H5Aopen(h5file, "timestep_fault", H5P_DEFAULT);
	checkH5Err(h5attr);
	checkH5Err(H5Aread(h5attr, H5T_NATIVE_INT, &timestepFault));
	checkH5Err(H5Aclose(h5attr));

	// Set the memory space (this is the same for all variables)
	hsize_t count[2] = {numSides(), numBndGP()};
	hid_t h5memSpace = H5Screate_simple(2, count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));

	// Offset for the file space
	hsize_t fStart[2] = {fileOffset(), 0};

	double* data[NUM_VARIABLES] = {mu, slipRate1, slipRate2, slip, slip1, slip2, state, strength};

	// Read the data
	for (unsigned int i = 0; i < NUM_VARIABLES; i++) {
		hid_t h5data = H5Dopen(h5file, VAR_NAMES[i], H5P_DEFAULT);
		checkH5Err(h5data);
		hid_t h5fSpace = H5Dget_space(h5data);
		checkH5Err(h5fSpace);

		// Read the data
		checkH5Err(H5Sselect_hyperslab(h5fSpace, H5S_SELECT_SET, fStart, 0L, count, 0L));

		checkH5Err(H5Dread(h5data, H5T_NATIVE_DOUBLE, h5memSpace, h5fSpace,
				h5XferList(), data[i]));

		checkH5Err(H5Sclose(h5fSpace));
		checkH5Err(H5Dclose(h5data));
	}

	checkH5Err(H5Sclose(h5memSpace));
	checkH5Err(H5Fclose(h5file));
}

void seissol::checkpoint::h5::Fault::write(int timestepFault)
{
	EPIK_TRACER("CheckPointFault_write");
	SCOREP_USER_REGION("CheckPointFault_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	if (numSides() == 0)
		return;

	logInfo(rank()) << "Checkpoint backend: Writing fault.";

	// Create array with all pointers
	EPIK_USER_REG(r_write_fault, "checkpoint_write_fault");
	SCOREP_USER_REGION_DEFINE(r_write_fault);
	EPIK_USER_START(r_write_fault);
	SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_fault", SCOREP_USER_REGION_TYPE_COMMON);

	// Attributes
	checkH5Err(H5Awrite(m_h5timestepFault[odd()], H5T_NATIVE_INT, &timestepFault));

	// Set memory and file space
	hsize_t fStart[2] = {fileOffset(), 0};
	hsize_t count[2] = {numSides(), numBndGP()};
	hid_t h5memSpace = H5Screate_simple(2, count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));
	checkH5Err(H5Sselect_hyperslab(m_h5fSpaceData, H5S_SELECT_SET, fStart, 0L, count, 0L));

	for (unsigned int i = 0; i < NUM_VARIABLES; i++) {
		checkH5Err(H5Dwrite(m_h5data[odd()][i], H5T_NATIVE_DOUBLE, h5memSpace, m_h5fSpaceData,
				h5XferList(), data(i)));
	}

	checkH5Err(H5Sclose(h5memSpace));

	EPIK_USER_END(r_write_fault);
	SCOREP_USER_REGION_END(r_write_fault);

	// Finalize the checkpoint
	finalizeCheckpoint();

	logInfo(rank()) << "Checkpoint backend: Writing fault. Done.";
}

bool seissol::checkpoint::h5::Fault::validate(hid_t h5file) const
{
	// Turn of error printing
	H5ErrHandler errHandler;

	// Check dimensions
	for (unsigned int i = 0; i < NUM_VARIABLES; i++) {
		hid_t h5data = H5Dopen(h5file, VAR_NAMES[i], H5P_DEFAULT);
		if (h5data < 0) {
			logWarning(rank()) << "Dataset" << VAR_NAMES[i] << "not found in checkpoint.";
			return false;
		}

		hid_t h5space = H5Dget_space(h5data);
		checkH5Err(H5Dclose(h5data));
		if (h5space < 0) {
			logWarning(rank()) << "Could not get space identifier for" << VAR_NAMES[i] << "in checkpoint.";
			return false;
		}

		bool isValid = true;

		int dims = H5Sget_simple_extent_ndims(h5space);
		if (dims != 2) {
			isValid = false;
			logWarning() << "Number of dimensions for" << VAR_NAMES[i] << "in checkpoint does not match.";
		} else {
			hsize_t dimSize[2];
			if (H5Sget_simple_extent_dims(h5space, dimSize, 0L) != 2) {
				isValid = false;
				logWarning(rank()) << "Could not get dimension sizes for" << VAR_NAMES[i] << "of checkpoint.";
			} else {
				if (dimSize[0] != numTotalElems()) {
					isValid = false;
					logWarning(rank()) << "Number of elements for" << VAR_NAMES[i] << "in checkpoint does not match.";
				}
				if (dimSize[1] != numBndGP()) {
					logWarning(rank()) << "Number of boundary points for" << VAR_NAMES[i] << "in checkpoint does not match.";
					isValid = false;
				}
			}
		}

		checkH5Err(H5Sclose(h5space));

		// Dimensions for at least one variable do not match -> invalid
		if (!isValid)
			return false;
	}

	return true;
}

hid_t seissol::checkpoint::h5::Fault::initFile(int odd, const char* filename)
{
	hid_t h5file;

	if (loaded()) {
		// Open the file
		h5file = open(filename, false);
		checkH5Err(h5file);

		// Fault writer
		m_h5timestepFault[odd] = H5Aopen(h5file, "timestep_fault", H5P_DEFAULT);
		checkH5Err(m_h5timestepFault[odd]);

		// Data
		for (unsigned int i = 0; i < NUM_VARIABLES; i++) {
			m_h5data[odd][i] = H5Dopen(h5file, VAR_NAMES[i], H5P_DEFAULT);
			checkH5Err(m_h5data[odd][i]);
		}
	} else {
		// Create the file
		hid_t h5plist = H5Pcreate(H5P_FILE_ACCESS);
		checkH5Err(h5plist);
		checkH5Err(H5Pset_libver_bounds(h5plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#ifdef USE_MPI
		MPIInfo info;
		checkH5Err(H5Pset_fapl_mpio(h5plist, comm(), info.get()));
#endif // USE_MPI

		h5file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, h5plist);
		checkH5Err(h5file);
		checkH5Err(H5Pclose(h5plist));

		// Create scalar dataspace for attributes
		hid_t h5spaceScalar = H5Screate(H5S_SCALAR);
		checkH5Err(h5spaceScalar);

		// Fault writer
		m_h5timestepFault[odd] = H5Acreate(h5file, "timestep_fault",
				H5T_STD_I32LE, h5spaceScalar, H5P_DEFAULT, H5P_DEFAULT);
		checkH5Err(m_h5timestepFault[odd]);
		int t = 0;
		checkH5Err(H5Awrite(m_h5timestepFault[odd], H5T_NATIVE_INT, &t));

		checkH5Err(H5Sclose(h5spaceScalar));

		// Variables
		for (unsigned int i = 0; i < NUM_VARIABLES; i++) {
			h5plist = H5Pcreate(H5P_DATASET_CREATE);
			checkH5Err(h5plist);
			checkH5Err(H5Pset_layout(h5plist, H5D_CONTIGUOUS));
			checkH5Err(H5Pset_alloc_time(h5plist, H5D_ALLOC_TIME_EARLY));
			m_h5data[odd][i] = H5Dcreate(h5file, VAR_NAMES[i], H5T_IEEE_F64LE, m_h5fSpaceData,
				H5P_DEFAULT, h5plist, H5P_DEFAULT);
			checkH5Err(m_h5data[odd][i]);
			checkH5Err(H5Pclose(h5plist));
		}
	}

	return h5file;
}


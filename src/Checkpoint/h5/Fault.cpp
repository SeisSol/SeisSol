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

#include "Fault.h"

#include <cassert>

#include "utils/stringutils.h"

void seissol::checkpoint::h5::Fault::setFilename(const char* filename)
{
	std::string file(filename);
	if (utils::StringUtils::endsWith(file, ".h5"))
		utils::StringUtils::replaceLast(file, ".h5", "-fault.h5");
	else
		file += "-fault.h5";

	CheckPoint::setFilename(file.c_str());
}

bool seissol::checkpoint::h5::Fault::init(
		double* mu, double* slipRate1, double* slipRate2, double* slip,
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
	m_mu = mu;
	m_slipRate1 = slipRate1;
	m_slipRate2 = slipRate2;
	m_slip = slip;
	m_state = state;
	m_strength = strength;
	m_numSides = numSides;
	m_numBndGP = numBndGP;

	// Compute total number of cells and local offset
	setSumOffset(numSides);

	// Dataspace for the file
	hsize_t fileSize[2] = {numTotalElems(), m_numBndGP};
	m_h5fSpaceData = H5Screate_simple(2, fileSize, 0L);

	setupXferList();

	return exists();
}

void seissol::checkpoint::h5::Fault::load(int &timestepFault)
{
	if (m_numSides == 0)
		return;

	logInfo(rank()) << "Loading fault checkpoint";

	hid_t h5file = open();
	checkH5Err(h5file);

	// Attributes
	hid_t h5attr = H5Aopen(h5file, "timestep_fault", H5P_DEFAULT);
	checkH5Err(h5attr);
	checkH5Err(H5Aread(h5attr, H5T_NATIVE_INT, &timestepFault));
	checkH5Err(H5Aclose(h5attr));

	// Create array with all pointers
	double * const data[] = {m_mu, m_slipRate1, m_slipRate2, m_slip, m_state, m_strength};
	assert(utils::ArrayUtils::size(data) == utils::ArrayUtils::size(VAR_NAMES));

	// Set the memory space (this is the same for all variables)
	hsize_t count[2] = {m_numSides, m_numBndGP};
	hid_t h5memSpace = H5Screate_simple(2, count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));

	// Offset for the file space
	hsize_t fStart[2] = {fileOffset(), 0};

	// Read the data
	for (unsigned int i = 0; i < utils::ArrayUtils::size(VAR_NAMES); i++) {
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

	if (m_numSides == 0)
		return;

	logInfo(rank()) << "Writing fault check point.";

	// Create array with all pointers
	const double * const data[] = {m_mu, m_slipRate1, m_slipRate2, m_slip, m_state, m_strength};
	assert(utils::ArrayUtils::size(data) == utils::ArrayUtils::size(VAR_NAMES));

	EPIK_USER_REG(r_write_fault, "checkpoint_write_fault");
	SCOREP_USER_REGION_DEFINE(r_write_fault);
	EPIK_USER_START(r_write_fault);
	SCOREP_USER_REGION_BEGIN(r_write_fault, "checkpoint_write_fault", SCOREP_USER_REGION_TYPE_COMMON);

	// Attributes
	checkH5Err(H5Awrite(m_h5timestepFault[odd()], H5T_NATIVE_INT, &timestepFault));

	// Set memory and file space
	hsize_t fStart[2] = {fileOffset(), 0};
	hsize_t count[2] = {m_numSides, m_numBndGP};
	hid_t h5memSpace = H5Screate_simple(2, count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));
	checkH5Err(H5Sselect_hyperslab(m_h5fSpaceData, H5S_SELECT_SET, fStart, 0L, count, 0L));

	for (unsigned int i = 0; i < utils::ArrayUtils::size(VAR_NAMES); i++) {
		checkH5Err(H5Dwrite(m_h5data[odd()][i], H5T_NATIVE_DOUBLE, h5memSpace, m_h5fSpaceData,
				h5XferList(), data[i]));
	}

	checkH5Err(H5Sclose(h5memSpace));

	EPIK_USER_END(r_write_fault);
	SCOREP_USER_REGION_END(r_write_fault);

	// Finalize the checkpoint
	finalizeCheckpoint();

	logInfo(rank()) << "Writing fault check point. Done.";
}

bool seissol::checkpoint::h5::Fault::validate(hid_t h5file) const
{
	// Turn of error printing
	H5ErrHandler errHandler;

	// Check dimensions
	for (unsigned int i = 0; i < utils::ArrayUtils::size(VAR_NAMES); i++) {
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
				if (dimSize[1] != m_numBndGP) {
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

hid_t seissol::checkpoint::h5::Fault::create(int odd, const char* filename)
{
	// Create the file
	hid_t h5plist = H5Pcreate(H5P_FILE_ACCESS);
	checkH5Err(h5plist);
	checkH5Err(H5Pset_libver_bounds(h5plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#ifdef USE_MPI
	checkH5Err(H5Pset_fapl_mpio(h5plist, comm(), MPI_INFO_NULL));
#endif // USE_MPI

	hid_t h5file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, h5plist);
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
	assert(utils::ArrayUtils::size(m_h5data[odd]) == utils::ArrayUtils::size(VAR_NAMES));

	for (unsigned int i = 0; i < utils::ArrayUtils::size(VAR_NAMES); i++) {
		h5plist = H5Pcreate(H5P_DATASET_CREATE);
		checkH5Err(h5plist);
		checkH5Err(H5Pset_layout(h5plist, H5D_CONTIGUOUS));
		checkH5Err(H5Pset_alloc_time(h5plist, H5D_ALLOC_TIME_EARLY));
		m_h5data[odd][i] = H5Dcreate(h5file, VAR_NAMES[i], H5T_IEEE_F64LE, m_h5fSpaceData,
			H5P_DEFAULT, h5plist, H5P_DEFAULT);
		checkH5Err(m_h5data[odd][i]);
		checkH5Err(H5Pclose(h5plist));
	}

	return h5file;
}

const char* seissol::checkpoint::h5::Fault::VAR_NAMES[6] = {
		"mu", "slip_rate1", "slip_rate2", "slip", "state", "strength" };


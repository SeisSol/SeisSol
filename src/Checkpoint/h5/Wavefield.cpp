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

#include "Parallel/MPI.h"

#include <cassert>

#include "utils/env.h"
#include "utils/mathutils.h"
#include "utils/stringutils.h"

#include "Wavefield.h"

#ifdef USE_MPI
#include "Checkpoint/MPIInfo.h"
#endif // USE_MPI

bool seissol::checkpoint::h5::Wavefield::init(size_t headerSize, unsigned long numDofs, unsigned int groupSize)
{
	seissol::checkpoint::Wavefield::init(headerSize, numDofs, groupSize);

	// Opaque data type for header
	// We cannot use compound data type on I/O nodes
	m_h5headerType = H5Tcreate(H5T_OPAQUE, headerSize);
	checkH5Err(m_h5headerType);

	// Data space for the file
	hsize_t fileSize = numTotalElems();
	m_h5fSpaceData = H5Screate_simple(1, &fileSize, 0L);
	checkH5Err(m_h5fSpaceData);

	setupXferList();

	return exists();
}

void seissol::checkpoint::h5::Wavefield::load(real* dofs)
{
	logInfo(rank()) << "Loading wave field checkpoint";

	seissol::checkpoint::CheckPoint::setLoaded();

	hid_t h5file = open(linkFile());
	checkH5Err(h5file);

	// Attributes
	hid_t h5attr = H5Aopen(h5file, "header", H5P_DEFAULT);
	checkH5Err(h5attr);
	checkH5Err(H5Aread(h5attr, m_h5headerType, header().data()));
	checkH5Err(H5Aclose(h5attr));

	// Get dataset
	hid_t h5data = H5Dopen(h5file, "values", H5P_DEFAULT);
	checkH5Err(h5data);
	hid_t h5fSpace = H5Dget_space(h5data);
	checkH5Err(h5fSpace);

	// Read the data
	unsigned int offset = 0;
	hsize_t fStart = fileOffset();
	hsize_t count = dofsPerIteration();
	hid_t h5memSpace = H5Screate_simple(1, &count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));
	for (unsigned int i = 0; i < totalIterations()-1; i++) {
		checkH5Err(H5Sselect_hyperslab(h5fSpace, H5S_SELECT_SET, &fStart, 0L, &count, 0L));

		checkH5Err(H5Dread(h5data, H5T_NATIVE_DOUBLE, h5memSpace, h5fSpace,
				h5XferList(), &dofs[offset]));

		// We are finished in less iterations, read data twice
		// so everybody needs the same number of iterations
		// TODO This does only work if everybody has at least 2 chunks
		if (i < iterations()-1) {
			fStart += count;
			offset += count;
		}
	}
	checkH5Err(H5Sclose(h5memSpace));

	// Read reminding data in the last iteration
	count = numDofs() - (iterations() - 1) * count;
	h5memSpace = H5Screate_simple(1, &count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));
	checkH5Err(H5Sselect_hyperslab(h5fSpace, H5S_SELECT_SET, &fStart, 0L, &count, 0L));
	checkH5Err(H5Dread(h5data, H5T_NATIVE_DOUBLE, h5memSpace, h5fSpace,
			h5XferList(), &dofs[offset]));
	checkH5Err(H5Sclose(h5memSpace));

	checkH5Err(H5Sclose(h5fSpace));
	checkH5Err(H5Dclose(h5data));
	checkH5Err(H5Fclose(h5file));
}

void seissol::checkpoint::h5::Wavefield::write(const void* header, size_t headerSize)
{
	EPIK_TRACER("CheckPoint_write");
	SCOREP_USER_REGION("CheckPoint_write", SCOREP_USER_REGION_TYPE_FUNCTION);

	logInfo(rank()) << "Checkpoint backend: Writing.";

	EPIK_USER_REG(r_header, "checkpoint_write_header");
	SCOREP_USER_REGION_DEFINE(r_header);
	EPIK_USER_START(r_header);
	SCOREP_USER_REGION_BEGIN(r_header, "checkpoint_write_header", SCOREP_USER_REGION_TYPE_COMMON);

	// Header
	checkH5Err(H5Awrite(m_h5header[odd()], m_h5headerType, header));

	EPIK_USER_END(r_header);
	SCOREP_USER_REGION_END(r_header);

	// Save data
	EPIK_USER_REG(r_write_wavefield, "checkpoint_write_wavefield");
	SCOREP_USER_REGION_DEFINE(r_write_wavefield);
	EPIK_USER_START(r_write_wavefield);
	SCOREP_USER_REGION_BEGIN(r_write_wavefield, "checkpoint_write_wavefield", SCOREP_USER_REGION_TYPE_COMMON);

	// Write the wave field
	unsigned int offset = 0;
	hsize_t fStart = fileOffset();
	hsize_t count = dofsPerIteration();
	hid_t h5memSpace = H5Screate_simple(1, &count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));
	for (unsigned int i = 0; i < totalIterations()-1; i++) {
		checkH5Err(H5Sselect_hyperslab(m_h5fSpaceData, H5S_SELECT_SET, &fStart, 0L, &count, 0L));

		checkH5Err(H5Dwrite(m_h5data[odd()], H5T_NATIVE_DOUBLE, h5memSpace, m_h5fSpaceData,
				h5XferList(), &const_cast<real*>(dofs())[offset]));

		// We are finished in less iterations, read data twice
		// so everybody needs the same number of iterations
		if (i < iterations()-1) {
			fStart += count;
			offset += count;
		}
	}
	checkH5Err(H5Sclose(h5memSpace));

	// Save reminding data in the last iteration
	count = numDofs() - (iterations() - 1) * count;
	h5memSpace = H5Screate_simple(1, &count, 0L);
	checkH5Err(h5memSpace);
	checkH5Err(H5Sselect_all(h5memSpace));
	checkH5Err(H5Sselect_hyperslab(m_h5fSpaceData, H5S_SELECT_SET, &fStart, 0L, &count, 0L));
	checkH5Err(H5Dwrite(m_h5data[odd()], H5T_NATIVE_DOUBLE, h5memSpace, m_h5fSpaceData,
			h5XferList(), &dofs()[offset]));
	checkH5Err(H5Sclose(h5memSpace));

	EPIK_USER_END(r_write_wavefield);
	SCOREP_USER_REGION_END(r_write_wavefield);

	// Finalize the checkpoint
	finalizeCheckpoint();

	logInfo(rank()) << "Checkpoint backend: Writing. Done.";
}

bool seissol::checkpoint::h5::Wavefield::validate(hid_t h5file) const
{
	// Turn of error printing
	H5ErrHandler errHandler;

	// Check #partitions
	hid_t h5attr = H5Aopen(h5file, "partitions", H5P_DEFAULT);
	if (h5attr < 0) {
		logWarning(rank()) << "Checkpoint does not have a partition attribute.";
		return false;
	}

	int p;
	herr_t err = H5Aread(h5attr, H5T_NATIVE_INT, &p);
	checkH5Err(H5Aclose(h5attr));
	if (err < 0 || p != partitions()) {
		logWarning(rank()) << "Partitions in checkpoint do not match.";
		return false;
	}

	// Check dimensions
	hid_t h5data = H5Dopen(h5file, "values", H5P_DEFAULT);
	if (h5data < 0) {
		logWarning(rank()) << "Checkpoint does not contains a data array.";
		return false;
	}

	hid_t h5space = H5Dget_space(h5data);
	checkH5Err(H5Dclose(h5data));
	if (h5space < 0) {
		logWarning(rank()) << "Could not get space identifier in checkpoint.";
		return false;
	}

	bool isValid = true;

	int dims = H5Sget_simple_extent_ndims(h5space);
	if (dims != 1) {
		isValid = false;
		logWarning(rank()) << "Number of dimensions in checkpoint does not match.";
	} else {
		hsize_t dimSize;
		if (H5Sget_simple_extent_dims(h5space, &dimSize, 0L) != 1) {
			isValid = false;
			logWarning(rank()) << "Could not get dimension sizes of checkpoint.";
		} else {
			if (dimSize != numTotalElems()) {
				isValid = false;
				logWarning(rank()) << "Number of elements in checkpoint does not match.";
			}
		}
	}
	checkH5Err(H5Sclose(h5space));

	return isValid;
}

hid_t seissol::checkpoint::h5::Wavefield::initFile(int odd, const char* filename)
{
	hid_t h5file;

	if (loaded()) {
		// Open the old file
		h5file = open(filename, false);
		checkH5Err(h5file);

		// Header
		m_h5header[odd] = H5Aopen(h5file, "header", H5P_DEFAULT);
		checkH5Err(m_h5header[odd]);

		// Data
		m_h5data[odd] = H5Dopen(h5file, "values", H5P_DEFAULT);
		checkH5Err(m_h5data[odd]);
	} else {
		// Create the file
		hid_t h5plist = H5Pcreate(H5P_FILE_ACCESS);
		checkH5Err(h5plist);
		checkH5Err(H5Pset_libver_bounds(h5plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
		hsize_t align = utils::Env::get<hsize_t>("SEISSOL_CHECKPOINT_ALIGNMENT", 0);
		if (align > 0)
			checkH5Err(H5Pset_alignment(h5plist, 1, align));
#ifdef USE_MPI
		MPIInfo info;
		checkH5Err(H5Pset_fapl_mpio(h5plist, seissol::MPI::mpi.comm(), info.get()));
#endif // USE_MPI

		h5file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, h5plist);
		checkH5Err(h5file);
		checkH5Err(H5Pclose(h5plist));

		// Create scalar dataspace for attributes
		hid_t h5spaceScalar = H5Screate(H5S_SCALAR);
		checkH5Err(h5spaceScalar);

		// Header
		m_h5header[odd] = H5Acreate(h5file, "header", m_h5headerType, h5spaceScalar,
				H5P_DEFAULT, H5P_DEFAULT);
		checkH5Err(m_h5header[odd]);

		// Partitions
		hid_t h5partitions = H5Acreate(h5file, "partitions", H5T_STD_I32LE, h5spaceScalar,
				H5P_DEFAULT, H5P_DEFAULT);
		checkH5Err(h5partitions);
		int p = partitions();
		checkH5Err(H5Awrite(h5partitions, H5T_NATIVE_INT, &p));
		checkH5Err(H5Aclose(h5partitions));

		checkH5Err(H5Sclose(h5spaceScalar));

		// Variable
		h5plist = H5Pcreate(H5P_DATASET_CREATE);
		checkH5Err(h5plist);
		checkH5Err(H5Pset_layout(h5plist, H5D_CONTIGUOUS));
		checkH5Err(H5Pset_alloc_time(h5plist, H5D_ALLOC_TIME_EARLY));
		m_h5data[odd] = H5Dcreate(h5file, "values", H5T_IEEE_F64LE, m_h5fSpaceData,
				H5P_DEFAULT, h5plist, H5P_DEFAULT);
		checkH5Err(m_h5data[odd]);
		checkH5Err(H5Pclose(h5plist));
	}

	return h5file;
}

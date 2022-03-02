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
 * Main class for MPI-IO checkpoints
 */

#ifndef CHECKPOINT_MPIO_CHECK_POINT_H
#define CHECKPOINT_MPIO_CHECK_POINT_H

#include <mpi.h>

#include <cassert>

#include "utils/env.h"

#include "Checkpoint/CheckPoint.h"
#include "Checkpoint/MPIInfo.h"
#include "Initializer/preProcessorMacros.fpp"
#include "SeisSol.h"

namespace seissol
{

namespace checkpoint
{

namespace mpio
{

class CheckPoint : virtual public seissol::checkpoint::CheckPoint
{
private:
	/** Identifiers of the files */
	MPI_File m_mpiFiles[2];

	/** True of the checkpoint files are opened */
	bool m_open;

	/** Number of bytes reserved for the header */
	unsigned long m_headerSize;

	/** The MPI data type for the header */
	MPI_Datatype m_headerType;

	/** The MPI data type for the file header */
	MPI_Datatype m_fileHeaderType;

	/** The MPI data type of the file data */
	MPI_Datatype m_fileDataType;

public:
	CheckPoint(unsigned long identifier)
		: seissol::checkpoint::CheckPoint(identifier),
		  m_open(false),
		  m_headerSize(0), m_headerType(MPI_DATATYPE_NULL),
		  m_fileHeaderType(MPI_DATATYPE_NULL), m_fileDataType(MPI_DATATYPE_NULL)
	{
	}

	virtual ~CheckPoint()
	{ }

	void setFilename(const char* filename)
	{
		initFilename(filename, "scp");
	}

	void close()
	{
		if (m_open) {
			for (unsigned int i = 0; i < 2; i++)
				MPI_File_close(&m_mpiFiles[i]);
			m_open = false;
		}

		if (m_headerType != MPI_DATATYPE_NULL) {
			MPI_Type_free(&m_headerType);
			MPI_Type_free(&m_fileDataType);
			if (rank() == 0)
				MPI_Type_free(&m_fileHeaderType);
		}
	}

protected:
	/**
	 * Submit the header data type
	 *
	 * @param type The MPI data type
	 */
	void setHeaderType(MPI_Datatype type)
	{
		m_headerType = type;
		MPI_Type_commit(&m_headerType);
	}

	/**
	 * Create the file view
	 *
	 * @param headerSize The size of the header in bytes
	 * @param elemSize The element size in bytes
	 * @param numElem The number of elements in the local part
	 */
	void defineFileView(unsigned long headerSize, unsigned int elemSize, unsigned long numElem, unsigned int numVars = 1)
	{
		// Check header size
		MPI_Aint lb, size;
		MPI_Type_get_extent(m_headerType, &lb, &size);
		if (size != static_cast<MPI_Aint>(headerSize))
			logError() << "Size of C struct and MPI data type do not match.";

		unsigned long align = utils::Env::get<unsigned long>("SEISSOL_CHECKPOINT_ALIGNMENT", 0);
		if (align > 0) {
			unsigned int blocks = (headerSize + align - 1) / align;
			m_headerSize = blocks * align;
		} else
			m_headerSize = headerSize;

		// Create element type
		MPI_Datatype elemType;
		MPI_Type_contiguous(elemSize, MPI_BYTE, &elemType);

		// Compute the number of blocks, we need to create the data type
		const unsigned long MAX_INT = 1ul<<30;
		unsigned int blocks = (numElem + MAX_INT - 1) / MAX_INT;

		// Create data file type
		int* blockLength = new int[blocks * numVars]; //{static_cast<int>(numElem)};
		MPI_Aint* displ = new MPI_Aint[blocks * numVars];
		unsigned long offset = m_headerSize + (fileOffset()-groupOffset()) * numVars * elemSize;
		for (unsigned int i = 0; i < numVars; i++) {
			// Jump to this process
			offset += groupOffset() * elemSize;

			for (unsigned int j = 0; j < blocks-1; j++) {
				blockLength[i*blocks + j] = MAX_INT;
				displ[i*blocks + j] = offset;
				offset += MAX_INT * elemSize;
			}

			blockLength[(i+1) * blocks - 1] = numElem - (blocks-1) * MAX_INT; // Set correct size for the last block
			displ[(i+1) * blocks - 1] = offset;
			offset += blockLength[(i+1) * blocks - 1] * elemSize;

			// Skip other processes after this in the group
			offset += (numGroupElems() - numElem - groupOffset()) * elemSize;
		}
		MPI_Type_create_hindexed(blocks * numVars, blockLength, displ, elemType, &m_fileDataType);
		MPI_Type_commit(&m_fileDataType);

		MPI_Type_free(&elemType);

		delete [] blockLength;
		delete [] displ;

		// Create the header file type
		if (rank() == 0) {
			MPI_Type_contiguous(headerSize, MPI_BYTE, &m_fileHeaderType);

			MPI_Type_commit(&m_fileHeaderType);
		} else
			// Only first rank write the header
			m_fileHeaderType = m_fileDataType;
	}

	bool exists()
	{
		if (!seissol::checkpoint::CheckPoint::exists())
			return false;

		MPI_File file = open();
		if (file == MPI_FILE_NULL)
			return false;

		bool hasCheckpoint = validate(file);
		MPI_File_close(&file);

		return hasCheckpoint;
	}

	void createFiles()
	{
		seissol::checkpoint::CheckPoint::createFiles();

		MPIInfo info;

		for (unsigned int i = 0; i < 2; i++) {
			checkMPIErr(MPI_File_open(comm(), const_cast<char*>(dataFile(i).c_str()),
					MPI_MODE_WRONLY | MPI_MODE_CREATE, info.get(), &m_mpiFiles[i]));

			// Sync file (required for performance measure)
			// TODO preallocate file first
			checkMPIErr(MPI_File_sync(m_mpiFiles[i]));
		}

		m_open = true;
	}

	/**
	 * Finalize checkpoint writing:
	 * Flush the file, update symbolic link, ...
	 */
	void finalizeCheckpoint()
	{
		EPIK_USER_REG(r_flush, "checkpoint_flush");
		SCOREP_USER_REGION_DEFINE(r_flush);
		EPIK_USER_START(r_flush);
		SCOREP_USER_REGION_BEGIN(r_flush, "checkpoint_flush", SCOREP_USER_REGION_TYPE_COMMON);

		checkMPIErr(MPI_File_sync(m_mpiFiles[odd()]));

		EPIK_USER_END(r_flush);
		SCOREP_USER_REGION_END(r_flush);
	}

	/**
	 * Open a check point file
	 *
	 * @return The MPI file handle
	 */
	MPI_File open()
	{
		MPI_File fh;

		MPIInfo info;

		int result = MPI_File_open(comm(), const_cast<char*>(linkFile()),
				MPI_MODE_RDONLY, info.get(), &fh);
		if (result != 0) {
			logWarning() << "Could not open checkpoint file";
			return MPI_FILE_NULL;
		}

		return fh;
	}

	/**
	 * Set header file view
	 *
	 * @return The MPI error code
	 */
	int setHeaderView(MPI_File file)
	{
		return MPI_File_set_view(file, 0, MPI_BYTE, m_fileHeaderType, const_cast<char*>("native"), MPI_INFO_NULL);
	}

	/**
	 * Set data file view
	 *
	 * @return The MPI error code
	 */
	int setDataView(MPI_File file)
	{
		return MPI_File_set_view(file, 0, MPI_BYTE, m_fileDataType, const_cast<char*>("native"), MPI_INFO_NULL);
	}

	/**
	 * @return The current MPI file
	 */
	MPI_File file() const
	{
		return m_mpiFiles[odd()];
	}

	/**
	 * @return The size of the header in bytes
	 */
	unsigned long headerSize() const
	{
		return m_headerSize;
	}

	MPI_Datatype headerType() const
	{
		return m_headerType;
	}

	/**
	 * Validate an existing check point file
	 */
	virtual bool validate(MPI_File file) = 0;

protected:
	static void checkMPIErr(int ret)
	{
		if (ret != 0) {
			char errString[MPI_MAX_ERROR_STRING+1];
			int length;
			MPI_Error_string(ret, errString, &length);
			assert(length < MPI_MAX_ERROR_STRING+1);

			errString[length] = '\0';
			logError() << "Error in the MPI checkpoint module:" << errString;
		}
	}
};

}

}

}

#endif // CHECKPOINT_MPIO_CHECK_POINT_H

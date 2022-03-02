/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Gilbert Brietzke (gilbert.brietzke AT lrz.de, http://www.lrz.de)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2015-2016, SeisSol Group
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
 * Main class for SIONlib checkpoints
 */

#ifndef CHECKPOINT_SIONLIB_CHECK_POINT_H
#define CHECKPOINT_SIONLIB_CHECK_POINT_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <errno.h>
#include <string>
#include <sys/stat.h>

#include <sion.h>

#include "utils/logger.h"
#include "utils/env.h"

#include "Checkpoint/CheckPoint.h"
#include "Initializer/preProcessorMacros.fpp"
#include "Initializer/typedefs.hpp"

namespace seissol
{

namespace checkpoint
{

namespace sionlib
{

class CheckPoint : virtual public seissol::checkpoint::CheckPoint
{
private:
	/** Number of files that should be used */
	int m_numFiles;

	/** (File system) Block size, -1 for auto */
	const sion_int32 m_blockSize;

	/** File name of the link */
	std::string m_linkFile;

	/** File name of odd/even checkpoint */
	std::string m_dataFile[2];

	/** Read mode */
	std::string m_readMode;

	/** Write mode */
	std::string m_writeMode;

	/** Rank of MPI_COMM_WORLD */
	int m_processIdentifier;

	/** The maximum written number of elements in chunk */
	unsigned int m_chunkElementCount;

public:
	CheckPoint(unsigned long identifier)
		: seissol::checkpoint::CheckPoint(identifier),
		m_numFiles(utils::Env::get<int>("SEISSOL_CHECKPOINT_SION_NUM_FILES", 1)),
		m_blockSize(utils::Env::get<sion_int32>("SEISSOL_CHECKPOINT_BLOCK_SIZE", -1)),
		m_readMode("rb"), m_writeMode("wb"),
		m_chunkElementCount(0)
	{
		// Detect ANSI/POSIX
		std::string backend = utils::Env::get<const char*>("SEISSOL_CHECKPOINT_SION_BACKEND", "ansi");
		utils::StringUtils::toLower(backend);
		if (backend == "ansi" || backend == "posix") {
			m_readMode += ","+backend;
			m_writeMode += ","+backend;
		}

		// Detect collective mode
		int collSize = utils::Env::get<int>("SEISSOL_CHECKPOINT_SION_COLL_SIZE", 0);
		if (collSize != 0) {
			std::string sCollSize = utils::StringUtils::toString(collSize);
			m_readMode += ",collective,collsize=" + sCollSize;
			m_writeMode += ",collective,collsize=" + sCollSize;

			std::string collMode = utils::Env::get<const char*>("SEISSOL_CHECKPOINT_SION_COLL_MODE", "merge");
			utils::StringUtils::toLower(collMode);
			if (collMode == "merge") {
				m_readMode += ",collectivemerge";
				m_writeMode += ",collectivemerge";
			}
		}

		MPI_Comm_rank(MPI_COMM_WORLD, &m_processIdentifier);
	}

	virtual ~CheckPoint()
	{}

	void setFilename(const char* filename)
	{
		initFilename(filename, 0L);

		m_linkFile = std::string(seissol::checkpoint::CheckPoint::linkFile())
		  + "/" + fname() + ".scp";
		for (unsigned int i = 0; i < 2; i++)
			m_dataFile[i] = seissol::checkpoint::CheckPoint::dataFile(i)
				+ "/" + fname() + ".scp";
	}

	/**
	 * Does nothing (all files are already closed)
	 */
	void close()
	{
	}

protected:
	void setChunkElementCount(unsigned int count)
	{
		m_chunkElementCount = count;
	}

	bool exists()
	{
		// Set the number of sion lib files
		// Needs to be done after initialization of the communicator but
		// before opening any files
		if (m_numFiles > partitions()) {
			logWarning(rank()) << "Reducing the number of SIONlib files for checkpointing to" << partitions();
			m_numFiles = partitions();
		}

		if (!seissol::checkpoint::CheckPoint::exists())
			return false;

		int file = open(linkFile(), m_readMode);
		if (file < 0) {
			logWarning() << "Could not open checkpoint file";
			return false;
		}

		int hasCheckpoint = validate(file);
		sionClose(file);

#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &hasCheckpoint, 1, MPI_INT, MPI_LAND, comm());
#endif // USE_MPI

		return hasCheckpoint;
	}

	void createFiles()
	{
		seissol::checkpoint::CheckPoint::createFiles();

		// Create the folder
		if (rank() == 0) {
			for (int i = 0; i < 2; i++) {
				int ret = mkdir(seissol::checkpoint::CheckPoint::dataFile(i).c_str(),
						S_IRWXU|S_IRWXG|S_IRWXO);
				if (ret < 0 && errno != EEXIST)
					checkErr(ret);
			}
		}

#ifdef USE_MPI // Make sure all processes see the folders
		MPI_Barrier(comm());
#endif // USE_MPI
	}

	/**
	 * @return The read mode used in SIONlib
	 */
	const std::string& readMode() const
	{
		return m_readMode;
	}

	/**
	 * @return The write mode used in SIONlib
	 */
	const std::string& writeMode() const
	{
		return m_writeMode;
	}

	/**
	 * Open a SIONlib file
	 *
	 * @param mode SIONlib open mode
	 * @return The file handle
	 */
	int open(const std::string& file, const std::string& mode)
	{
		int numFiles = m_numFiles;
		MPI_Comm localComm = comm();
		sion_int64 chunksize = m_chunkElementCount * sizeof(real);
		sion_int32 fsblksize = m_blockSize;
		int gRank = m_processIdentifier;
		return sion_paropen_mpi(file.c_str(), mode.c_str(), &numFiles, comm(), &localComm,
				&chunksize, &fsblksize, &gRank, 0L, 0L);
	}

	void finalizeCheckpoint(int file)
	{
		SCOREP_USER_REGION_DEFINE(r_flush);
		SCOREP_USER_REGION_BEGIN(r_flush, "checkpoint_flush", SCOREP_USER_REGION_TYPE_COMMON);

		sionClose(file);

		SCOREP_USER_REGION_END(r_flush);
	}

	const std::string linkFile() const
	{
		return m_linkFile;
	}

	const std::string dataFile(int odd) const
	{
		return m_dataFile[odd];
	}

	/**
	 * Close a SIONlib file
	 *
	 * @param file The file handle
	 */
	void sionClose(int file)
	{
		checkErr(sion_parclose_mpi(file), SION_SUCCESS);
	}

private:
	/** Validate an existing check point file */
	bool validate(int file) const
	{
		unsigned long id;
		size_t size = sion_coll_fread(&id, sizeof(id), 1, file);
		if (size != 1) {
			logWarning() << "Could not read checkpoint header";
			return false;
		}

		if (id != identifier()) {
			logWarning() << "Checkpoint identifier does not match";
			return false;
		}

		return true;
	}

protected:
	template<typename T>
	static void checkErr(T ret) {
		if (ret < 0)
			logError() << "Error in the SIONlib checkpoint module";
	}

	/** Can be used to check read/write errors
	 * @param ret The return value
	 * @param target The expected return value (> 0)
	 */
	template<typename T, typename U>
	static void checkErr(T ret, U target) {
		checkErr(ret);
		if (ret != target)
			logError() << "Error in the SIONlib checkpoint module: Expected:"
				<< target << "Received:" << ret;
	}
};

}

}

}

#endif // CHECKPOINT_SIONLIB_CHECK_POINT_H

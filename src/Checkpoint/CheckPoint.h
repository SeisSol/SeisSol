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
 * Common interface for checkpoints
 */

#ifndef CHECKPOINT_CHECK_POINT_H
#define CHECKPOINT_CHECK_POINT_H

#ifdef USE_MPI
#include <mpi.h>
#endif // USE_MPI

#include <cstdio>
#include <string>
#include <unistd.h>
#include <sys/stat.h>

#include "utils/logger.h"
#include "utils/path.h"
#include "utils/stringutils.h"

namespace seissol
{

namespace checkpoint
{

/**
 * Common interface for all checkpoints
 */
class CheckPoint
{
private:
	/** Checkpoint identifier (written to the beginning of the file) */
	const unsigned long m_identifier;

	/** File name of the symbolic link to the current checkpoint */
	std::string m_linkFile;

	/** Path to the HDF5 files */
	utils::Path m_dir;

	/** File name of the HDF5 files (without path) */
	std::string m_files[2];

	/** The rank of the process */
	int m_rank;

	/** Total number of partitions */
	int m_partitions;

#ifdef USE_MPI
	/** Communicator used for this checkpoint */
	MPI_Comm m_comm;
#endif // USE_MPI

	/** Next checkpoint should go to even or odd? */
	int m_odd;

	/** Total number of elements */
	unsigned long m_numTotalElems;

	/** Offset in the file (in elements) */
	unsigned long m_fileOffset;

	/** Number of ranks per group */
	unsigned int m_groupSize;

	/** Total number of elements in the group */
	unsigned long m_numGroupElems;

	/** Offset in the group */
	unsigned long m_groupOffset;

	/** Was the checkpoint loaded */
	bool m_loaded;

public:
	CheckPoint(unsigned long identifier)
		: m_identifier(identifier),
		  m_rank(0), m_partitions(1), // default for no MPI
#ifdef USE_MPI
		  m_comm(MPI_COMM_NULL),
#endif // USE_MPI
		  m_odd(0), // Start with even checkpoint
		  m_numTotalElems(0), m_fileOffset(0),
		  m_groupSize(0), m_numGroupElems(0), m_groupOffset(0),
		  m_loaded(false)
	{}

	virtual ~CheckPoint() {}

	/**
	 * Set the file names for the output files and the symbolic link
	 *
	 * @param filename Base file name
	 */
	virtual void setFilename(const char* filename) = 0;

	/**
	 * Should be called when a checkpoint is loaded
	 */
	void setLoaded()
	{
		m_loaded = true;
	}

	/**
	 * Update checkpoint symlink
	 */
	virtual void updateLink()
	{
		// Update the link to the latest checkpoint file
		if (m_rank == 0) {
			remove(linkFile());
			if (symlink(m_files[m_odd].c_str(), linkFile()) != 0)
				logWarning() << "Failed to create symbolic link to current checkpoint.";
		}

		// Switch between even and odd
		m_odd = 1 - m_odd;
	}

	/**
	 * Close the checkpoint files
	 */
	virtual void close() = 0;

protected:
	/**
	 * Initializes file names and link file name
	 *
	 * @param filename The base file name
	 * @param extension The extension
	 */
	virtual void initFilename(const char* filename, const char* extension)
	{
		std::string ext;
		if (extension)
			ext = std::string(".") + extension;

		// Link file
		m_linkFile = filename;
		if (extension && !utils::StringUtils::endsWith(m_linkFile, ext))
			m_linkFile += ext;

		m_dir = utils::Path(m_linkFile).dir();

		// Real output files
		m_files[0] = m_files[1] = utils::Path(m_linkFile).basename();
		if (extension) {
			std::string ext0 = ".0" + ext;
			utils::StringUtils::replaceLast(m_files[0], ext, ext0);
			std::string ext1 = ".1" + ext;
			utils::StringUtils::replaceLast(m_files[1], ext, ext1);
		} else {
			m_files[0] += ".0";
			m_files[1] += ".1";
		}
	}

#ifdef USE_MPI
	/**
	 * Sets up the MPI communicator
	 *
	 * @param comm The MPI communicator
	 */
	void setComm(MPI_Comm comm)
	{
		m_comm = comm;

		MPI_Comm_rank(m_comm, &m_rank);
		MPI_Comm_size(m_comm, &m_partitions);
	}
#endif // USE_MPI

	/**
	 * Compute the total number of elements and the local offset in the checkpoint file
	 *
	 * @param numElems Local number of elements
	 */
	void setSumOffset(unsigned long numElems)
	{
		m_numTotalElems = numElems;
		m_fileOffset = numElems;
#ifdef USE_MPI
		MPI_Allreduce(MPI_IN_PLACE, &m_numTotalElems, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm());
		MPI_Scan(MPI_IN_PLACE, &m_fileOffset, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm());
#endif // USE_MPI
		m_fileOffset -= numElems;
	}

	void setGroupSumOffset(unsigned long numElems, unsigned int groupSize = 1)
	{
		m_groupSize = groupSize;
		m_numGroupElems = numElems;
		m_groupOffset = numElems;

		if (groupSize > 1) {
#ifdef USE_MPI
			// Create a communicator for the group to do the reductions
			MPI_Comm groupComm;
			MPI_Comm_split(comm(), m_rank / groupSize, m_rank, &groupComm);

			MPI_Allreduce(MPI_IN_PLACE, &m_numGroupElems, 1, MPI_UNSIGNED_LONG, MPI_SUM, groupComm);
			MPI_Scan(MPI_IN_PLACE, &m_groupOffset, 1, MPI_UNSIGNED_LONG, MPI_SUM, groupComm);

			MPI_Comm_free(&groupComm);
#endif // USE_MPI
		}

		m_groupOffset -= numElems;
	}

	/**
	 * Checks if a valid checkpoint exists
	 *
	 * @return True if a checkpoint exists, false otherwise
	 */
	bool exists() const
	{
		struct stat statBuffer;
		if (lstat(m_linkFile.c_str(), &statBuffer) != 0)
			// No link exists
			return false;

		if (stat(m_linkFile.c_str(), &statBuffer) != 0) {
			logWarning(m_rank) << "Link exists, but checkpoint file was not found.";
			return false;
		}

		return true;
	}

	/**
	 * Creates the new checkpoint files. Override this in the backend
	 * @warning This might delete any existing checkpoint from previous executions
	 */
	virtual void createFiles()
	{
		// If a checkpoint was loaded, we will reopen the old files
		// otherwise, create a backup
		if (!m_loaded) {
			for (unsigned int i = 0; i < 2; i++) {
				std::string filename = dataFile(i);
				// Backup checkpoint
				if (m_rank == 0)
					rename(filename.c_str(), (filename + ".bak").c_str());
			}

#ifdef USE_MPI
			// Make sure the file is moved before anyone create a new file
			MPI_Barrier(m_comm);
#endif // USE_MPI
		}
	}

	/**
	 * @return File name of the symbolic link
	 */
	const char* linkFile() const
	{
		return m_linkFile.c_str();
	}

	/**
	 * @return Filename of the odd or even data file
	 */
	std::string dataFile(int odd) const
	{
		return m_dir + utils::Path(m_files[odd]);
	}

	/**
	 * @return The name for files if the back-end uses directories
	 */
	virtual const char* fname() const
	{
		return "cp";
	}

	/**
	 * @return The identifier of the file
	 */
	unsigned long identifier() const
	{
		return m_identifier;
	}

	int rank() const
	{
		return m_rank;
	}

	int partitions() const
	{
		return m_partitions;
	}

	bool loaded() const
	{
		return m_loaded;
	}

#ifdef USE_MPI
	MPI_Comm comm() const
	{
		return m_comm;
	}
#endif // USE_MPI

	int odd() const
	{
		return m_odd;
	}

	unsigned long numTotalElems() const
	{
		return m_numTotalElems;
	}

	unsigned long fileOffset() const
	{
		return m_fileOffset;
	}

	unsigned int groupSize() const
	{
		return m_groupSize;
	}

	unsigned long numGroupElems() const
	{
		return m_numGroupElems;
	}

	unsigned long groupOffset() const
	{
		return m_groupOffset;
	}
};

}

}

#endif // CHECKPOINT_CHECK_POINT_H

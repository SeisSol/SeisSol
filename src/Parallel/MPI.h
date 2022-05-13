/**
 * @file
 * This file is part of SeisSol.
 *
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
 * MPI Wrapper
 */

#ifndef MPI_H
#define MPI_H

#ifndef USE_MPI
#include "MPIDummy.h"
#else // USE_MPI

#include <mpi.h>

#include "utils/logger.h"

#include "MPIBasic.h"

#ifdef ACL_DEVICE
#include <cstdlib>
#include <string>
#include <sstream>
#include <device.h>
#endif  // ACL_DEVICE

#endif // USE_MPI

namespace seissol
{

#ifndef USE_MPI
typedef MPIDummy MPI;
#else // USE_MPI

/**
 * MPI handling.
 *
 * Make sure only one instance of this class exists!
 */
class MPI : public MPIBasic
{
private:
	MPI_Comm m_comm;

#ifdef ACL_DEVICE
    int m_localRank{};
    int m_localSize{};
    int m_deviceId{};
#endif // ACL_DEVICE

private:
	MPI()
		: m_comm(MPI_COMM_NULL)
	{ }

public:
	~MPI()
	{ }

#ifdef ACL_DEVICE

  private:
    /**
     * @brief Reads and returns environment variables
     *
     * Some MPI vendors usually provides env. variables which allows to find out the local rank and size
     * before calling MPI_Init(...). However, they tend to name these variables differently, i.e. uniquely
     * for their implementation. Thus, the function take some potential candidates and loop through them and try
     * to retrieve a value.
     *
     * @param candidates a vector of strings with names of possible env. variables
     * @throws std::string in case if a value cannot get retrieved from a candidate list
     * @throws std::invalid_argument in case if an env. variable doesn't contain an integer, e.g. char, string, etc.
     * @throws std::out_of_range in case is an env. variable contains a value bigger that a size of integer
     * */
    static int readValueFromEnvVariables(std::vector<std::string> &candidates) {
      char* valueStr = nullptr;
      for (auto envVar: candidates) {
        valueStr = std::getenv(envVar.c_str());
        if (valueStr)
          break;
      }

      if (!valueStr) {
        std::stringstream stream;

        stream << "could not detect any env. variable from a list of candidates, namely: ";
        for (const auto& item: candidates) {
          stream << item << ", ";
        }
        stream << ". Please, consider to use any other MPI implementation with an offloading support.";

        logError() << stream.str();
      }

      return std::stoi(std::string(valueStr));
    }

  public:
    /**
     * @brief Inits Device(s).
     *
     * Some MPI implementations create a so-called context between GPUs and OS Processes inside of MPI_Init(...).
     * It results in allocating some memory buffers in memory attached to the nearest NUMA domain
     * of a core where a process is running. In case of somebody wants to bind a processes in a different way,
     * e.g. move a process closer to a GPU, it must be done before calling MPI_Init(...) using env. variables
     * or hwloc library.
     *
     * Currently, the function does a simple binding, i.e. a local rank controls the corresponding devices.
     * For instance, localRank=2 is going to use deviceId=2. The user is responsible for the correct binding.
     * She/he must refer to a documentation of their job scheduler or MPI implementation to achieve correct
     * GPU/CPU affinity. Note, one can improve the current binding strategy using hwloc.
     * See, Professional CUDA programming, subsection Affinity on MPI-CUDA Programs as a reference.
     *
     * The function supports the following MPI implementations: OpenMPI, MVAPICH2, IntelMPI
     * */
    void  bindRankToDevice() {
      try {
        std::vector<std::string> rankEnvVars{{"OMPI_COMM_WORLD_LOCAL_RANK"},
                                             {"MV2_COMM_WORLD_LOCAL_RANK"},
                                             {"SLURM_LOCALID"}, {"PMI_RANK"} };

        std::vector<std::string> sizeEnvVars{{"OMPI_COMM_WORLD_LOCAL_SIZE"},
                                             {"MV2_COMM_WORLD_LOCAL_SIZE"},
                                             {"SLURM_NTASKS_PER_NODE"}, {"PMI_SIZE"}};

        m_localRank = readValueFromEnvVariables(rankEnvVars);
        m_localSize = readValueFromEnvVariables(sizeEnvVars);
      }
      catch (const std::invalid_argument &err) {
        logError() << err.what() << ". File: " << __FILE__ << ", line: " << __LINE__;
      }
      catch (const std::out_of_range& err) {
        logError() << err.what() << ". File: " << __FILE__ << ", line: " << __LINE__;
      }

      device::DeviceInstance& device = device::DeviceInstance::getInstance();
      int m_numDevices = device.api->getNumDevices();
      if (m_localSize > m_numDevices) {
        logError() << "Local mpi size (in a compute node) is greater than the number of avaliable devices."
                   << "Over-subscription of devices is currently not supported in Seissol."
                   << "Adjust num. local mpi rank and num. local devices.\n"
                   << "File: " << __FILE__ << ", line: " << __LINE__;
      }
      m_deviceId = m_localRank;

#ifdef _OPENMP
#pragma omp parallel
      {
#pragma omp critical
        {
          device.api->setDevice(m_deviceId);
        }
      }
#else
      device.api->setDevice(m_deviceId);
#endif
    }
    int getDeviceID() { return m_deviceId; }
#endif // ACL_DEVICE

	/**
	 * Initialize MPI
	 */
	void init(int &argc, char** &argv)
	{
	  // Note: Strictly speaking, we only require MPI_THREAD_MULTIPLE if using
	  // a communication thread and/or async I/O.
	  // The safer (and more sane) option is to enable it by default.
		int required = MPI_THREAD_MULTIPLE;
		int provided;
		MPI_Init_thread(&argc, &argv, required, &provided);

		setComm(MPI_COMM_WORLD);

		// Test this after setComm() to get the correct m_rank
		if (required < provided)
			logWarning(m_rank) << utils::nospace << "Required MPI thread support (" << required
				<< ") is smaller than provided thread support (" << provided << ").";
	}

	void setComm(MPI_Comm comm)
	{
		m_comm = comm;

		MPI_Comm_rank(comm, &m_rank);
		MPI_Comm_size(comm, &m_size);
	}

	/**
	 * @return The main communicator for the application
	 */
	MPI_Comm comm() const
	{
		return m_comm;
	}

	void barrier(MPI_Comm comm) const
	{
		MPI_Barrier(comm);
	}

	/**
	 * Finalize MPI
	 */
	void finalize()
	{
		fault.finalize();

		MPI_Finalize();
	}

public:
	/** The only instance of the class */
	static MPI mpi;
};

#endif // USE_MPI

}

#endif // MPI_H

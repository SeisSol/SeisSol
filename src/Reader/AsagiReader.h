/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2016-2017, SeisSol Group
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
 * Velocity field reader Fortran interface
 */

#ifndef ASAGIREADER_H
#define ASAGIREADER_H

#include "Parallel/MPI.h"
#include <asagi.h>
#include <easi/util/AsagiReader.h>

#include "utils/env.h"
#include "utils/logger.h"

#include "AsagiModule.h"
#include "Monitoring/instrumentation.fpp"

namespace seissol
{
namespace asagi
{
enum NUMACache_Mode
{
	NUMA_OFF, NUMA_ON, NUMA_CACHE
};

class AsagiReader : public easi::AsagiReader
{
private:
	/** Prefix for environment variables */
	const std::string m_envPrefix;

	/** Number of threads used by ASAGI */
	unsigned int m_asagiThreads;
  
#ifdef USE_MPI
  /** MPI communicator used by ASAGI */
  MPI_Comm m_comm;
#endif

public:
	AsagiReader(  const char* envPrefix
#ifdef USE_MPI
                , MPI_Comm comm = seissol::MPI::mpi.comm()
#endif                
                ) : m_envPrefix(envPrefix)
#ifdef USE_MPI
                    , m_comm(comm)
#endif
  {}
  
  virtual ::asagi::Grid* open(char const* file, char const* varname);
  virtual unsigned numberOfThreads() const { return m_asagiThreads; }

private:
	static NUMACache_Mode getNUMAMode();
};

/**
 *
 * @param file File name of the netCDF file
 * @param varname The variable name in the netCDF file
 * @return ASAGI grid
 */
::asagi::Grid* AsagiReader::open(const char* file, const char* varname) {
  SCOREP_USER_REGION("AsagiReader_open", SCOREP_USER_REGION_TYPE_FUNCTION);

  const int rank = seissol::MPI::mpi.rank();

  ::asagi::Grid* grid = ::asagi::Grid::createArray();

  if (utils::Env::get<bool>((m_envPrefix + "_SPARSE").c_str(), false)) {
    grid->setParam("GRID", "CACHE");
  }

  // Set MPI mode
  if (AsagiModule::mpiMode() != MPI_OFF) {
#ifdef USE_MPI
    ::asagi::Grid::Error err = grid->setComm(m_comm);
    if (err != ::asagi::Grid::SUCCESS)
      logError() << "Could not set ASAGI communicator:" << err;

#endif // USE_MPI

    if (AsagiModule::mpiMode() == MPI_COMM_THREAD)
      grid->setParam("MPI_COMMUNICATION", "THREAD");
  }

  // Set NUMA mode
  m_asagiThreads = utils::Env::get((m_envPrefix  + "_NUM_THREADS").c_str(), 0u);
  if (m_asagiThreads == 0)
    m_asagiThreads = AsagiModule::totalThreads();
  else if (static_cast<int>(m_asagiThreads) > AsagiModule::totalThreads()) {
    logWarning(rank) << "Only" << AsagiModule::totalThreads()
        << "threads can be used for ASAGI initialization.";
    m_asagiThreads = AsagiModule::totalThreads();
  }

  if (AsagiModule::mpiMode() == MPI_COMM_THREAD)
    m_asagiThreads--; // one thread is used for communication

  grid->setThreads(m_asagiThreads);

  switch (getNUMAMode()) {
  case NUMA_ON:
    grid->setParam("NUMA_COMMUNICATION", "ON");
    break;
  case NUMA_OFF:
    grid->setParam("NUMA_COMMUNICATION", "OFF");
    break;
  case NUMA_CACHE:
    grid->setParam("NUMA_COMMUNICATION", "CACHE");
    break;
  }

  // Set vertex centered grid
  grid->setParam("VALUE_POSITION", "VERTEX_CENTERED");

  // Set additional parameters
  std::string blockSize = utils::Env::get((m_envPrefix  + "_BLOCK_SIZE").c_str(), "64");
  grid->setParam("BLOCK_SIZE_0", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_1", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_2", blockSize.c_str());

  std::string cacheSize = utils::Env::get((m_envPrefix  + "_CACHE_SIZE").c_str(), "128");
  grid->setParam("CACHE_SIZE", cacheSize.c_str());

  grid->setParam("VARIABLE", varname);

  bool abort = false;
  // Read the data
  //SCOREP_RECORDING_OFF();
#ifdef _OPENMP
  #pragma omp parallel shared(abort) num_threads(m_asagiThreads)
#endif // _OPENMP
  {
    ::asagi::Grid::Error err = grid->open(file);
    if (err != ::asagi::Grid::SUCCESS)
      abort = true;
  }
  //SCOREP_RECORDING_ON();
  if (abort) {
    logError() << "Could not open " << file << " with ASAGI.";
    return nullptr;
  }

  return grid;
}

NUMACache_Mode AsagiReader::getNUMAMode()
{
  const char* numaModeName = utils::Env::get("SEISSOL_ASAGI_NUMA_MODE", "ON");

  if (strcmp(numaModeName, "ON") == 0)
    return NUMA_ON;
  if (strcmp(numaModeName, "OFF") == 0)
    return NUMA_OFF;
  if (strcmp(numaModeName, "CACHE") == 0)
    return NUMA_CACHE;

  logError() << "Unknown NUMA mode:" << numaModeName;
  return NUMA_OFF;
}
}
}

#endif // ASAGIREADER_H

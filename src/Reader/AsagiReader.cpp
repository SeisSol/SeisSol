// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifdef USE_ASAGI

#include "AsagiReader.h"

#include <Monitoring/Instrumentation.h>
#include <Reader/AsagiModule.h>
#include <asagi.h>
#include <mpi.h>
#include <utils/env.h>
#include <utils/logger.h>

namespace seissol::asagi {
/**
 *
 * @param file File name of the netCDF file
 * @param varname The variable name in the netCDF file
 * @return ASAGI grid
 */
::asagi::Grid* AsagiReader::open(const char* file, const char* varname) {
  SCOREP_USER_REGION("AsagiReader_open", SCOREP_USER_REGION_TYPE_FUNCTION);

  ::asagi::Grid* grid = ::asagi::Grid::createArray();

  if (utils::Env::get<bool>(envPrefix + "_SPARSE", false)) {
    grid->setParam("GRID", "CACHE");
  }

  // Set MPI mode
  if (AsagiModule::mpiMode() != AsagiMPIMode::Off) {
#ifdef USE_MPI
    ::asagi::Grid::Error const err = grid->setComm(comm);
    if (err != ::asagi::Grid::SUCCESS) {
      logError() << "Could not set ASAGI communicator:" << err;
    }

#endif // USE_MPI

    if (AsagiModule::mpiMode() == AsagiMPIMode::CommThread) {
      grid->setParam("MPI_COMMUNICATION", "THREAD");
    }
  }

  // Set NUMA mode
  asagiThreads = utils::Env::get(envPrefix + "_NUM_THREADS", 0U);
  if (asagiThreads == 0) {
    asagiThreads = AsagiModule::totalThreads();
  } else if (static_cast<int>(asagiThreads) > AsagiModule::totalThreads()) {
    logWarning() << "Only" << AsagiModule::totalThreads()
                 << "threads can be used for ASAGI initialization.";
    asagiThreads = AsagiModule::totalThreads();
  }

  if (AsagiModule::mpiMode() == AsagiMPIMode::CommThread) {
    // one thread is used for communication
    --asagiThreads;
  }

  grid->setThreads(asagiThreads);

  switch (getNumaMode()) {
  case NumaCacheMode::On:
    grid->setParam("NUMA_COMMUNICATION", "ON");
    break;
  case NumaCacheMode::Off:
    grid->setParam("NUMA_COMMUNICATION", "OFF");
    break;
  case NumaCacheMode::Cache:
    grid->setParam("NUMA_COMMUNICATION", "CACHE");
    break;
  }

  // Set vertex centered grid
  grid->setParam("VALUE_POSITION", "VERTEX_CENTERED");

  // Set additional parameters
  const std::string blockSize = utils::Env::get(envPrefix + "_BLOCK_SIZE", "64");
  grid->setParam("BLOCK_SIZE_0", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_1", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_2", blockSize.c_str());

  const std::string cacheSize = utils::Env::get(envPrefix + "_CACHE_SIZE", "128");
  grid->setParam("CACHE_SIZE", cacheSize.c_str());

  grid->setParam("VARIABLE", varname);

  bool abort = false;
  // Read the data
  // SCOREP_RECORDING_OFF();
#ifdef _OPENMP
#pragma omp parallel shared(abort) num_threads(asagiThreads)
#endif // _OPENMP
  {
    const ::asagi::Grid::Error err = grid->open(file);
    if (err != ::asagi::Grid::SUCCESS) {
      abort = true;
    }
  }
  // SCOREP_RECORDING_ON();
  if (abort) {
    logError() << "Could not open " << file << " with ASAGI.";
    return nullptr;
  }

  return grid;
}

NumaCacheMode AsagiReader::getNumaMode() {
  const std::string numaModeName = utils::Env::get("SEISSOL_ASAGI_NUMA_MODE", "ON");

  if (numaModeName == "ON") {
    return NumaCacheMode::On;
  }
  if (numaModeName == "OFF") {
    return NumaCacheMode::Off;
  }
  if (numaModeName == "CACHE") {
    return NumaCacheMode::Cache;
  }

  logError() << "Unknown NUMA mode:" << numaModeName;
  return NumaCacheMode::Off;
}

unsigned AsagiReader::numberOfThreads() const { return asagiThreads; }

AsagiReader::AsagiReader(const char* envPrefix
#ifdef USE_MPI
                         ,
                         MPI_Comm comm
#endif
                         )
    : envPrefix(envPrefix)
#ifdef USE_MPI
      ,
      comm(comm)
#endif
{
}

} // namespace seissol::asagi

#endif

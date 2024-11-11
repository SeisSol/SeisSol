#ifdef USE_ASAGI

#include "AsagiReader.h"

namespace seissol::asagi {
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

  if (utils::Env::get<bool>((envPrefix + "_SPARSE").c_str(), false)) {
    grid->setParam("GRID", "CACHE");
  }

  // Set MPI mode
  if (AsagiModule::mpiMode() != AsagiMPIMode::Off) {
#ifdef USE_MPI
    ::asagi::Grid::Error err = grid->setComm(comm);
    if (err != ::asagi::Grid::SUCCESS)
      logError() << "Could not set ASAGI communicator:" << err;

#endif // USE_MPI

    if (AsagiModule::mpiMode() == AsagiMPIMode::CommThread)
      grid->setParam("MPI_COMMUNICATION", "THREAD");
  }

  // Set NUMA mode
  asagiThreads = utils::Env::get((envPrefix + "_NUM_THREADS").c_str(), 0u);
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
  std::string blockSize = utils::Env::get((envPrefix + "_BLOCK_SIZE").c_str(), "64");
  grid->setParam("BLOCK_SIZE_0", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_1", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_2", blockSize.c_str());

  std::string cacheSize = utils::Env::get((envPrefix + "_CACHE_SIZE").c_str(), "128");
  grid->setParam("CACHE_SIZE", cacheSize.c_str());

  grid->setParam("VARIABLE", varname);

  bool abort = false;
  // Read the data
  // SCOREP_RECORDING_OFF();
#ifdef _OPENMP
#pragma omp parallel shared(abort) num_threads(asagiThreads)
#endif // _OPENMP
  {
    ::asagi::Grid::Error err = grid->open(file);
    if (err != ::asagi::Grid::SUCCESS)
      abort = true;
  }
  // SCOREP_RECORDING_ON();
  if (abort) {
    logError() << "Could not open " << file << " with ASAGI.";
    return nullptr;
  }

  return grid;
}

NumaCacheMode AsagiReader::getNumaMode() {
  const char* numaModeName = utils::Env::get("SEISSOL_ASAGI_NUMA_MODE", "ON");

  if (strcmp(numaModeName, "ON") == 0)
    return NumaCacheMode::On;
  if (strcmp(numaModeName, "OFF") == 0)
    return NumaCacheMode::Off;
  if (strcmp(numaModeName, "CACHE") == 0)
    return NumaCacheMode::Cache;

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

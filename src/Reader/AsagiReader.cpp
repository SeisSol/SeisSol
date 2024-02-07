#include "AsagiReader.h"

#ifdef USE_ASAGI

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
  m_asagiThreads = utils::Env::get((m_envPrefix + "_NUM_THREADS").c_str(), 0u);
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
  std::string blockSize = utils::Env::get((m_envPrefix + "_BLOCK_SIZE").c_str(), "64");
  grid->setParam("BLOCK_SIZE_0", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_1", blockSize.c_str());
  grid->setParam("BLOCK_SIZE_2", blockSize.c_str());

  std::string cacheSize = utils::Env::get((m_envPrefix + "_CACHE_SIZE").c_str(), "128");
  grid->setParam("CACHE_SIZE", cacheSize.c_str());

  grid->setParam("VARIABLE", varname);

  bool abort = false;
  // Read the data
  // SCOREP_RECORDING_OFF();
#ifdef _OPENMP
#pragma omp parallel shared(abort) num_threads(m_asagiThreads)
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

NUMACache_Mode AsagiReader::getNUMAMode() {
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

unsigned AsagiReader::numberOfThreads() const { return m_asagiThreads; }

AsagiReader::AsagiReader(const char* envPrefix
#ifdef USE_MPI
                         ,
                         MPI_Comm comm
#endif
                         )
    : m_envPrefix(envPrefix)
#ifdef USE_MPI
      ,
      m_comm(comm)
#endif
{
}

} // namespace seissol::asagi

#endif

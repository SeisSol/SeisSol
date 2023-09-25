#include "ResultWriter/ThreadsPinningWriter.h"
#include "Parallel/MPI.h"
#include "Common/filesystem.h"
#include <sched.h>
#include <fstream>

#ifndef __APPLE__
#include <sys/sysinfo.h>
#ifdef USE_NUMA_AWARE_PINNING
#include <numa.h>
#endif // USE_NUMA_AWARE_PINNING
#endif // __APPLE__

#ifndef __APPLE__
namespace seissol::writer::pinning::details {
struct PinningInfo {
  std::string coreIds{};
  std::string numaIds{};
};

PinningInfo getPinningInfo(cpu_set_t const& set) {
  std::stringstream coreIdsStream;
  std::stringstream numaIdsStream;
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    if (CPU_ISSET(cpu, &set)) {
      coreIdsStream << cpu << ',';
#ifdef USE_NUMA_AWARE_PINNING
      numaIdsStream << numa_node_of_cpu(cpu) << ',';
#endif // USE_NUMA_AWARE_PINNING
    }
  }

  auto trim = [](std::string& str) {
    if (str.empty()) {
#ifdef USE_NUMA_AWARE_PINNING
      str = std::string("no-info");
#else
      str = std::string("none");
#endif // USE_NUMA_AWARE_PINNING
    } else {
      str.pop_back();
    }
  };

  PinningInfo pinningInfo;
  pinningInfo.coreIds = coreIdsStream.str();
  pinningInfo.numaIds = numaIdsStream.str();

  trim(pinningInfo.coreIds);
  trim(pinningInfo.numaIds);

  return pinningInfo;
}
} // namespace seissol::writer::pinning::details
#endif // __APPLE__

void seissol::writer::ThreadsPinningWriter::write(const seissol::parallel::Pinning& pinning) {
  logWarning(MPI::mpi.rank()) << "ThreadsPinningWriter is not supported on MacOS.";
#ifndef __APPLE__
  auto workerInfo = pinning::details::getPinningInfo(pinning.getWorkerUnionMask().set);
#ifdef USE_COMM_THREAD
  auto freeCpus = pinning.getFreeCPUsMask();
  auto commThreadInfo = pinning::details::getPinningInfo(freeCpus);
#else
  cpu_set_t emptyUnion;
  CPU_ZERO(&emptyUnion);
  auto commThreadInfo = pinning::details::getPinningInfo(emptyUnion);

#endif // USE_COMM_THREAD

  auto workerThreads = seissol::MPI::mpi.collectContainer(workerInfo.coreIds);
  auto workerNumas = seissol::MPI::mpi.collectContainer(workerInfo.numaIds);

  auto commThreads = seissol::MPI::mpi.collectContainer(commThreadInfo.coreIds);
  auto commNumas = seissol::MPI::mpi.collectContainer(commThreadInfo.numaIds);

  auto localRanks = seissol::MPI::mpi.collect(seissol::MPI::mpi.sharedMemMpiRank());
  auto numNProcs = seissol::MPI::mpi.collect(get_nprocs());

  if (seissol::MPI::mpi.rank() == 0) {
    seissol::filesystem::path path(outputDirectory);
    path += seissol::filesystem::path("-threadPinning.csv");

    std::fstream fileStream(path, std::ios::out);
    fileStream
        << "hostname,rank,localRank,workermask,workernuma,commthread_mask,commthread_numa,nproc\n";

    const auto& hostNames = seissol::MPI::mpi.getHostNames();
    for (int rank = 0; rank < seissol::MPI::mpi.size(); ++rank) {
      fileStream << "\"" << hostNames[rank] << "\"," << rank << ',' << localRanks[rank] << ",\""
                 << workerThreads[rank] << "\",\"" << workerNumas[rank] << "\",\""
                 << commThreads[rank] << "\",\"" << commNumas[rank] << "\"," << numNProcs[rank]
                 << "\n";
    }

    fileStream.close();
  }
#endif // __APPLE__
}

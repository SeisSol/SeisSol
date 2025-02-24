// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ResultWriter/ThreadsPinningWriter.h"
#include "Common/Filesystem.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include <Parallel/Pin.h>
#include <fstream>
#include <ios>
#include <sched.h>
#include <sstream>
#include <string>

#ifndef __APPLE__
#include <sys/sysinfo.h>
#ifdef USE_NUMA_AWARE_PINNING
#include <numa.h>
#endif // USE_NUMA_AWARE_PINNING
#endif // __APPLE__

#ifndef __APPLE__
namespace {

using namespace seissol::parallel;

struct PinningInfo {
  std::string coreIds;
  std::string numaIds;
};

PinningInfo getPinningInfo(const cpu_set_t& set) {
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
} // namespace
#endif // __APPLE__

void seissol::writer::ThreadsPinningWriter::write(const seissol::parallel::Pinning& pinning) {
#ifndef __APPLE__
  auto workerInfo = getPinningInfo(seissol::parallel::Pinning::getWorkerUnionMask().set);

  PinningInfo commThreadInfo;
  if (seissol::useCommThread(seissol::MPI::mpi)) {
    auto freeCpus = pinning.getFreeCPUsMask();
    commThreadInfo = getPinningInfo(freeCpus.set);
  } else {
    cpu_set_t emptyUnion;
    CPU_ZERO(&emptyUnion);
    commThreadInfo = getPinningInfo(emptyUnion);
  }

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
#else
  logWarning() << "ThreadsPinningWriter is not supported on MacOS.";
#endif // __APPLE__
}

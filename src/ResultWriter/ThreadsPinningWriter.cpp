// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ResultWriter/ThreadsPinningWriter.h"

#include "Common/Filesystem.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <fstream>
#include <ios>
#include <sched.h>
#include <sstream>
#include <string>
#include <utils/env.h>

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

void seissol::writer::ThreadsPinningWriter::write(const seissol::parallel::Pinning& pinning,
                                                  utils::Env& env) {
#ifndef __APPLE__
  auto workerInfo = getPinningInfo(seissol::parallel::Pinning::getWorkerUnionMask().set);

  PinningInfo commThreadInfo;
  if (seissol::useCommThread(seissol::Mpi::mpi, env)) {
    auto freeCpus = pinning.getFreeCPUsMask();
    commThreadInfo = getPinningInfo(freeCpus.set);
  } else {
    cpu_set_t emptyUnion;
    CPU_ZERO(&emptyUnion);
    commThreadInfo = getPinningInfo(emptyUnion);
  }

  auto workerThreads = seissol::Mpi::mpi.collectContainer(workerInfo.coreIds);
  auto workerNumas = seissol::Mpi::mpi.collectContainer(workerInfo.numaIds);

  auto commThreads = seissol::Mpi::mpi.collectContainer(commThreadInfo.coreIds);
  auto commNumas = seissol::Mpi::mpi.collectContainer(commThreadInfo.numaIds);

  auto localRanks = seissol::Mpi::mpi.collect(seissol::Mpi::mpi.sharedMemMpiRank());
  auto numNProcs = seissol::Mpi::mpi.collect(get_nprocs());

  if (seissol::Mpi::mpi.rank() == 0) {
    std::filesystem::path path(outputDirectory);
    path += std::filesystem::path("-threadPinning.csv");

    std::fstream fileStream(path, std::ios::out);
    fileStream << "hostname,device,rank,localRank,workermask,workernuma,commthread_mask,commthread_"
                  "numa,nproc\n";

    const auto& hostNames = seissol::Mpi::mpi.getHostNames();
    const auto& pcis = seissol::Mpi::mpi.getPCIAddresses();
    const std::string nullstring;
    for (int rank = 0; rank < seissol::Mpi::mpi.size(); ++rank) {
      const auto& pci = pcis.empty() ? nullstring : pcis[rank];
      fileStream << "\"" << hostNames[rank] << "\",\"" << pci << "\"," << rank << ','
                 << localRanks[rank] << ",\"" << workerThreads[rank] << "\",\"" << workerNumas[rank]
                 << "\",\"" << commThreads[rank] << "\",\"" << commNumas[rank] << "\","
                 << numNProcs[rank] << "\n";
    }

    fileStream.close();
  }
#else
  logWarning() << "ThreadsPinningWriter is not supported on MacOS.";
#endif // __APPLE__
}

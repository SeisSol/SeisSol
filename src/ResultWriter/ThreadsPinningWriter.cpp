// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ResultWriter/ThreadsPinningWriter.h"

#include "Common/Filesystem.h"
#include "IO/Instance/Point/Csv.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <algorithm>
#include <cstddef>
#include <sched.h>
#include <sstream>
#include <string>
#include <utils/env.h>
#include <vector>

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
    const auto& hostNames = seissol::Mpi::mpi.getHostNames();
    const auto& pcis = seissol::Mpi::mpi.getPCIAddresses();
    const std::string nullstring;

    // the masks are a list of the cores a rank runs on, so how wide they get is a property of the
    // machine rather than something to guess at
    const auto widest = [](const std::vector<std::string>& values) {
      std::size_t width = 1;
      for (const auto& value : values) {
        width = std::max(width, value.size());
      }
      return width;
    };

    seissol::io::instance::point::Csv table("threadPinning");
    table.addTextColumn("hostname", widest(hostNames));
    table.addTextColumn("device", widest(pcis));
    table.addColumn<int>("rank");
    table.addColumn<int>("localRank");
    table.addTextColumn("workermask", widest(workerThreads));
    table.addTextColumn("workernuma", widest(workerNumas));
    table.addTextColumn("commthread_mask", widest(commThreads));
    table.addTextColumn("commthread_numa", widest(commNumas));
    table.addColumn<int>("nproc");

    for (int rank = 0; rank < seissol::Mpi::mpi.size(); ++rank) {
      table.addText(hostNames[rank]);
      table.addText(pcis.empty() ? nullstring : pcis[rank]);
      table.addCell<int>(rank);
      table.addCell<int>(localRanks[rank]);
      table.addText(workerThreads[rank]);
      table.addText(workerNumas[rank]);
      table.addText(commThreads[rank]);
      table.addText(commNumas[rank]);
      table.addCell<int>(numNProcs[rank]);
    }

    seissol::filesystem::path path(outputDirectory_);
    path += seissol::filesystem::path("-threadPinning.csv");
    table.writeFile(path.string());
  }
#else
  logWarning() << "ThreadsPinningWriter is not supported on MacOS.";
#endif // __APPLE__
}

// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Lukas Krenz

#include "Pin.h"

#include "Parallel/MPI.h"
#include "utils/logger.h"
#include <Common/IntegerMaskParser.h>
#include <async/as/Pin.h>
#include <cassert>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <mpi.h>
#include <sched.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#ifndef __APPLE__
#include <sys/sysinfo.h>
#ifdef USE_NUMA_AWARE_PINNING
#include "numa.h"
#endif // USE_NUMA_AWARE_PINNING
#endif // __APPLE__

namespace seissol::parallel {

using namespace async::as;

std::deque<bool> Pinning::parseOnlineCpuMask(std::string s, unsigned numberOfConfiguredCpus) {
  std::deque<bool> onlineMask(numberOfConfiguredCpus, false);

  // The file has the format e.g. 0-1,12-59
  // Possibly also on some systems 0,12
  // or just one range 0-7

  // Step 1: Split into tokens
  // E.g. 0-1, 12-59 into t1=0-1, t2=12-59
  std::vector<std::string> tokens;
  size_t pos = 0;
  std::string token;
  while ((pos = s.find(',')) != std::string::npos) {
    token = s.substr(0, pos);
    s.erase(0, pos + 1);
    tokens.push_back(token);
  }
  tokens.push_back(s);

  // Step 2: Set mask for each token
  for (auto& t : tokens) {
    pos = t.find('-');
    int beginRange = 0;
    int endRange = 0;
    if (pos == std::string::npos) {
      beginRange = std::stoi(t);
      endRange = beginRange;
    } else {
      beginRange = std::stoi(t.substr(0, pos));
      t.erase(0, pos + 1);
      endRange = std::stoi(t);
    }

    for (int cpu = beginRange; cpu <= endRange; ++cpu) {
      onlineMask[cpu] = true;
    }
  }
  return onlineMask;
}

CpuMask seissol::parallel::Pinning::computeOnlineCpuMask() {
#ifndef __APPLE__
  CPU_ZERO(&onlineMask.set);
  std::deque<bool> mask;

  const std::string onlineFilePath = "/sys/devices/system/cpu/online";
  const std::ifstream file(onlineFilePath);

  if (file.good()) {
    std::stringstream buffer;
    buffer << file.rdbuf();
    mask = parseOnlineCpuMask(buffer.str(), get_nprocs_conf());

  } else {
    logWarning() << "Could not read" << onlineFilePath << "Assuming that all cpus are online.";
    mask = std::deque<bool>(get_nprocs_conf(), true);
  }

  assert(static_cast<int>(mask.size()) == get_nprocs_conf());
  for (unsigned cpu = 0; cpu < mask.size(); ++cpu) {
    if (mask[cpu]) {
      CPU_SET(cpu, &onlineMask.set);
    }
  }
  return CpuMask{onlineMask};
#else
  return {};
#endif
}

Pinning::Pinning() {
  // Affinity mask for the OpenMP workers
  openmpMask = getWorkerUnionMask();
  computeOnlineCpuMask();
}

void Pinning::checkEnvVariables() {
#ifndef __APPLE__
  if (const char* envVariable = std::getenv("SEISSOL_FREE_CPUS_MASK")) {
    auto parsedResult = seissol::IntegerMaskParser::parse(std::string(envVariable));
    if (parsedResult) {
      parsedFreeCPUsMask = parsedResult.value();

      bool isMaskGood{true};
      const auto numLocalProcesses = MPI::mpi.sharedMemMpiSize();
      if (numLocalProcesses > static_cast<int>(parsedFreeCPUsMask.size())) {
        logInfo() << "There are more communication (and/or output-writing) threads"
                  << "to pin than locations defined in `SEISSOL_FREE_CPUS_MASK`";

        isMaskGood = false;
      } else {
        const auto maxCpuId = get_nprocs();
        for (auto localProcessId = 0; localProcessId < static_cast<int>(parsedFreeCPUsMask.size());
             ++localProcessId) {
          for (auto cpu : parsedFreeCPUsMask[localProcessId]) {
            if (cpu > maxCpuId) {
              logInfo() << "Free cpu mask of the local process" << localProcessId
                        << "is out of bounds. CPU/core id" << cpu << "exceeds max. value"
                        << maxCpuId;
              isMaskGood = false;
              break;
            }
          }
        }
      }

      if (isMaskGood) {
        logInfo() << "Binding free cpus according to `SEISSOL_FREE_CPUS_MASK` env. variable.";
      } else {
        logWarning() << "Ignoring `SEISSOL_FREE_CPUS_MASK` env. variable.";
        logWarning() << "`SEISSOL_FREE_CPUS_MASK` Format:"
                     << "(<int>|<range: int-int>|<list: {int,+}>),+";
        parsedFreeCPUsMask = IntegerMaskParser::MaskType{};
      }
    } else {
      logWarning() << "Failed to parse `SEISSOL_FREE_CPUS_MASK` env. variable";
    }
  }
#endif // __APPLE__
}

CpuMask Pinning::getWorkerUnionMask() {
#ifndef __APPLE__
  cpu_set_t workerUnion;
  CPU_ZERO(&workerUnion);
#ifdef _OPENMP
#pragma omp parallel default(none) shared(workerUnion)
  {
    cpu_set_t worker;
    CPU_ZERO(&worker);
    sched_getaffinity(0, sizeof(cpu_set_t), &worker);
#pragma omp critical
    {
      CPU_OR(&workerUnion, &workerUnion, &worker);
    }
  }
#else
  sched_getaffinity(0, sizeof(cpu_set_t), &workerUnion);
#endif

  return CpuMask{workerUnion};

#else
  return CpuMask{};
#endif // __APPLE__
}

CpuMask Pinning::getFreeCPUsMask() const {
#ifndef __APPLE__
  const auto nodeOpenMpMask = getNodeMask();

  cpu_set_t freeMask{};
  CPU_ZERO(&freeMask);

  if (not parsedFreeCPUsMask.empty()) {
    const auto localProcessor = MPI::mpi.sharedMemMpiRank();
    for (const auto& cpu : parsedFreeCPUsMask[localProcessor]) {
      CPU_SET(cpu, &freeMask);
    }
    return CpuMask{freeMask};
  }

#ifdef USE_NUMA_AWARE_PINNING
  // Find all numa nodes on which some OpenMP worker is pinned to
  std::set<int> numaDomainsOfThisProcess{};
  for (int cpu = 0; cpu < get_nprocs_conf(); ++cpu) {
    if (CPU_ISSET(cpu, &openmpMask.set)) {
      numaDomainsOfThisProcess.insert(numa_node_of_cpu(cpu));
    }
  }

  // Set free mask to all free threads which are on one of our numa nodes
  for (int cpu = 0; cpu < get_nprocs_conf(); ++cpu) {
    const bool isOnline = CPU_ISSET(cpu, &onlineMask.set);
    const bool isFree = !CPU_ISSET(cpu, &nodeOpenMpMask.set);
    if (isOnline && isFree) {
      const int numaNode = numa_node_of_cpu(cpu);
      const bool isValidNumaNode = numaDomainsOfThisProcess.count(numaNode) != 0;
      if (isValidNumaNode) {
        CPU_SET(cpu, &freeMask);
      }
    }
  }
#else
  // Set now contains all unused cores on the machine.
  // Note that pinning of the communication thread is then not Numa-aware if there's more than one
  // rank per node!
  for (int cpu = 0; cpu < get_nprocs_conf(); ++cpu) {
    const bool isOnline = CPU_ISSET(cpu, &onlineMask.set);
    const bool isFree = !CPU_ISSET(cpu, &nodeOpenMpMask.set);
    if (isOnline && isFree) {
      CPU_SET(cpu, &freeMask);
    }
  }
#endif // USE_NUMA_AWARE_PINNING

  return CpuMask{freeMask};
#else
  return {};
#endif // __APPLE__
}

bool Pinning::freeCPUsMaskEmpty(const CpuMask& mask) {
#ifndef __APPLE__
  return CPU_COUNT(&(mask.set)) == 0;
#else
  return false;
#endif // __APPLE__
}

CpuMask Pinning::getOnlineMask() const { return onlineMask; }

bool Pinning::areAllCpusOnline() {
#ifndef __APPLE__
  return get_nprocs_conf() == get_nprocs();
#else
  return true;
#endif
}

void Pinning::pinToFreeCPUs() const {
  auto freeMask = getFreeCPUsMask();
#ifndef __APPLE__
  sched_setaffinity(0, sizeof(cpu_set_t), &(freeMask.set));
#endif // __APPLE__
}

std::string Pinning::maskToString(const CpuMask& mask) {
#ifndef __APPLE__
  const auto& set = mask.set;
  std::stringstream st;
  for (int cpu = 0; cpu < get_nprocs_conf(); ++cpu) {
    if (cpu % 10 == 0 && cpu != 0 && cpu != get_nprocs_conf() - 1) {
      st << '|';
    }
    if (CPU_ISSET(cpu, &set)) {
      st << cpu % 10;
    } else {
      st << '-';
    }
  }
  return st.str();

#else
  return "Affinity is not supported on MacOS.";
#endif // __APPLE__
}

CpuMask Pinning::getNodeMask() {
#ifndef __APPLE__
  const auto workerMask = getWorkerUnionMask().set;

  // We have to use this due to the insanity of std::vector<bool>
  auto workerMaskArray = std::vector<char>(get_nprocs_conf(), 0);
  for (int cpu = 0; cpu < get_nprocs_conf(); ++cpu) {
    workerMaskArray[cpu] = CPU_ISSET(cpu, &workerMask);
  }

  MPI_Allreduce(MPI_IN_PLACE,
                workerMaskArray.data(),
                workerMaskArray.size(),
                MPI_CHAR,
                MPI_BOR,
                MPI::mpi.sharedMemComm());

  cpu_set_t nodeMask;
  CPU_ZERO(&nodeMask);
  for (int cpu = 0; cpu < get_nprocs_conf(); ++cpu) {
    const auto isSet = workerMaskArray[cpu] != 0;
    if (isSet) {
      CPU_SET(cpu, &nodeMask);
    }
  }

  return CpuMask{nodeMask};
#else
  return {};
#endif // __APPLE__
}

} // namespace seissol::parallel

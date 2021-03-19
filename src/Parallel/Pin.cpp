/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * 
 * @section LICENSE
 * Copyright (c) 2019, SeisSol Group
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
 * 
 **/

#include "Pin.h"

#include <sys/sysinfo.h>
#include <sched.h>
#include <sstream>
#include <set>
#include "Parallel/MPI.h"

#ifdef USE_NUMA_AWARE_PINNING
#include "numa.h"
#endif

seissol::parallel::Pinning::Pinning() {
  // Affinity mask for the OpenMP workers
  openmpMask = getWorkerUnionMask();
}

cpu_set_t seissol::parallel::Pinning::getWorkerUnionMask() const {
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

  return workerUnion;
}

cpu_set_t seissol::parallel::Pinning::getFreeCPUsMask() const {
  const auto nodeOpenMpMask = getNodeMask();

  cpu_set_t freeMask{};
  CPU_ZERO(&freeMask);

#ifdef USE_NUMA_AWARE_PINNING
  // Find all numa nodes on which some OpenMP worker is pinned to
  std::set<int> numaDomainsOfThisProcess{};
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    if (CPU_ISSET(cpu, &openmpMask)) {
      numaDomainsOfThisProcess.insert(numa_node_of_cpu(cpu));
    }
  }

  // Set free mask to all free threads which are on one of our numa nodes
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    const bool isFree = !CPU_ISSET(cpu, &nodeOpenMpMask);
    if (isFree) {
      const int numaNode = numa_node_of_cpu(cpu);
      const bool isValidNumaNode = numaDomainsOfThisProcess.count(numaNode) != 0;
      if (isValidNumaNode) {
        CPU_SET(cpu, &freeMask);
      }
    }
  }
#else
  // Set now contains all unused cores on the machine.
  // Note that pinning of the communication thread is then not Numa-aware if there's more than one rank per node!
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    if (!CPU_ISSET(cpu, &nodeOpenMpMask)) {
      CPU_SET(cpu, &freeMask);
    }
  }
#endif


  return freeMask;
}

bool seissol::parallel::Pinning::freeCPUsMaskEmpty(cpu_set_t const& set) {
  return CPU_COUNT(&set) == 0;
}

void seissol::parallel::Pinning::pinToFreeCPUs() const {
  auto freeMask = getFreeCPUsMask();
  sched_setaffinity(0, sizeof(cpu_set_t), &freeMask);
}

std::string seissol::parallel::Pinning::maskToString(cpu_set_t const& set) {
  std::stringstream st;
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    if (cpu % 10 == 0 && cpu != 0 && cpu != get_nprocs()-1) {
      st << '|';
    }
    if (CPU_ISSET(cpu, &set)) {
      st << cpu % 10;
    } else {
      st << '-';
    }
  }
  return st.str();
}

cpu_set_t seissol::parallel::Pinning::getNodeMask() const {
  const auto workerMask = getWorkerUnionMask();

  // We have to use this due to the insanity of std::vector<bool>
  auto workerMaskArray = std::vector<char>( get_nprocs(), 0);
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    workerMaskArray[cpu] = CPU_ISSET(cpu, &workerMask);
  }

  MPI_Comm commNode;
  MPI_Comm_split_type(MPI::mpi.comm(), MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &commNode);
  MPI_Allreduce(MPI_IN_PLACE, workerMaskArray.data(), workerMaskArray.size(), MPI_CHAR, MPI_BOR, commNode);

  cpu_set_t nodeMask;
  CPU_ZERO(&nodeMask);
  for (int cpu = 0; cpu < get_nprocs(); ++cpu) {
    const auto isSet = workerMaskArray[cpu] != 0;
    if (isSet) {
      CPU_SET(cpu, &nodeMask);
    }
  }

  return nodeMask;
}

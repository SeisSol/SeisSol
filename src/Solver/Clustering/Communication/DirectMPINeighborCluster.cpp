// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DirectMPINeighborCluster.h"
#include <Memory/MemoryAllocator.h>
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/MPI.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <memory>
#include <mpi.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

#include "utils/logger.h"

namespace seissol::solver::clustering::communication {

bool DirectMPINeighborClusterUnidirectional::poll() {
  if (requests.size() > 0) {
    std::lock_guard guard(requestMutex);
    if (enqueueCounterTrue != enqueueCounterPoll) {
      MPI_Startall(requests.size(), requests.data());
      ++enqueueCounterPoll;
    }
    if (enqueueCounterPoll != dequeueCounterPoll) {
      int result = 0;
      MPI_Testall(requests.size(), requests.data(), &result, MPI_STATUSES_IGNORE);
      const auto done = result != 0;
      if (done) {
        ++dequeueCounterPoll;
      }
    }
    return dequeueCounterPoll == enqueueCounter;
  } else {
    return true;
  }
}

void DirectMPINeighborClusterUnidirectional::start(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    {
      std::lock_guard guard(requestMutex);
      assert(dequeueCounterPoll == enqueueCounter);
      ++enqueueCounter;
    }
    runtime.enqueueHost([&] {
      std::lock_guard guard(requestMutex);
      ++enqueueCounterTrue;
    });
  }
}

DirectMPINeighborClusterUnidirectional::~DirectMPINeighborClusterUnidirectional() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
}

} // namespace seissol::solver::clustering::communication

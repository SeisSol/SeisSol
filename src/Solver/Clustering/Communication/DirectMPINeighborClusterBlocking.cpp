// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DirectMPINeighborClusterBlocking.h"
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

bool DirectMPINeighborClusterBlocking::poll() {
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

void DirectMPINeighborClusterBlocking::start(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    std::lock_guard guard(requestMutex);
    assert(dequeueCounterPoll == enqueueCounter);
    runtime.enqueueHost([&] {
      std::lock_guard guard(requestMutex);
      ++enqueueCounterTrue;
    });
    ++enqueueCounter;
  }
}
DirectMPISendNeighborClusterBlocking::DirectMPISendNeighborClusterBlocking(
    const std::vector<RemoteCluster>& remote,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority)
    : SendNeighborCluster(cpuExecutor, priority) {
  requests.resize(remote.size());
  for (std::size_t i = 0; i < remote.size(); ++i) {
    MPI_Send_init(remote[i].data,
                  remote[i].size,
                  remote[i].datatype,
                  remote[i].rank,
                  remote[i].tag,
                  MPI::mpi.comm(),
                  &requests[i]);
  }
}

DirectMPINeighborClusterBlocking::~DirectMPINeighborClusterBlocking() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
}

DirectMPIRecvNeighborClusterBlocking::DirectMPIRecvNeighborClusterBlocking(
  const std::vector<RemoteCluster>& remote,
  const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
  double priority)
  : RecvNeighborCluster(cpuExecutor, priority) {
requests.resize(remote.size());
for (std::size_t i = 0; i < remote.size(); ++i) {
  MPI_Recv_init(remote[i].data,
                remote[i].size,
                remote[i].datatype,
                remote[i].rank,
                remote[i].tag,
                MPI::mpi.comm(),
                &requests[i]);
}
}

bool DirectMPISendNeighborClusterBlocking::poll() {
  return DirectMPINeighborClusterBlocking::poll();
}

bool DirectMPIRecvNeighborClusterBlocking::poll() {
  return DirectMPINeighborClusterBlocking::poll();
}

void DirectMPISendNeighborClusterBlocking::start(parallel::runtime::StreamRuntime& runtime) {
  DirectMPINeighborClusterBlocking::start(runtime);
}

void DirectMPIRecvNeighborClusterBlocking::start(parallel::runtime::StreamRuntime& runtime) {
  DirectMPINeighborClusterBlocking::start(runtime);
}

void DirectMPISendNeighborClusterBlocking::stop(parallel::runtime::StreamRuntime& runtime) {
}

void DirectMPIRecvNeighborClusterBlocking::stop(parallel::runtime::StreamRuntime& runtime) {
}

/*
bool CopyMPISendNeighborCluster::poll() {
    return true;
}

void CopyMPISendNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
    void* splitEvent = runtime.recordEvent();
    for (void* stream : copyStreams) {
        device::DeviceInstance::getInstance().api->syncStreamWithEvent(stream, splitEvent);
        device::DeviceInstance::getInstance().api->copyToAsync(TODO, TODO, TODO, stream);
        device::DeviceInstance::getInstance().api->streamHostFunction(stream, []() {
            MPI_Start(TODO);
        });
        device::DeviceInstance::getInstance().api->streamHostFunction(stream, []() {
            MPI_Wait(TODO, MPI_STATUS_IGNORE);
        });
    }
    for (void* stream : copyStreams) {
        void* joinEvent = runtime.getEvent();
        device::DeviceInstance::getInstance().api->recordEventOnStream(joinEvent, stream);
        device::DeviceInstance::getInstance().api->syncStreamWithEvent(runtime.stream(), joinEvent);
    }
#endif


#ifdef ACL_DEVICE
    runtime.enqueueHost([=]() {
        MPI_Startall(TODO, TODO);
    });
    void* splitEvent = runtime.recordEvent();
    for (void* stream : copyStreams) {
        device::DeviceInstance::getInstance().api->syncStreamWithEvent(stream, splitEvent);
        device::DeviceInstance::getInstance().api->streamHostFunction(stream, []() {
            MPI_Wait(TODO, MPI_STATUS_IGNORE);
        });
        device::DeviceInstance::getInstance().api->copyToAsync(TODO, TODO, TODO, stream);
    }
    for (void* stream : copyStreams) {
        void* joinEvent = runtime.getEvent();
        device::DeviceInstance::getInstance().api->recordEventOnStream(joinEvent, stream);
        device::DeviceInstance::getInstance().api->syncStreamWithEvent(runtime.stream(), joinEvent);
    }
#endif
    // runtime.stream()
    // split into streams
    // TODO: copy here
    // then run direct MPI logic, only MPI_Start
    // reunite
}

bool CopyMPIRecvNeighborCluster::poll() {
    return true;
}

for CopyMPIRecv:
MPI_startall
split into streams
MPI_Wait per stream
copy per stream
reunite
*/

} // namespace seissol::solver::clustering::communication

// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DirectMPINeighborCluster.h"
#include <Parallel/MPI.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <mpi.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::solver::clustering::communication {

bool DirectMPISendNeighborCluster::poll() {
  if (requests.size() > 0) {
    std::lock_guard guard(requestMutex);
    int result;
    MPI_Testsome(requests.size(), requests.data(), &result, status.data(), MPI_STATUSES_IGNORE);
    return result == requests.size();
  } else {
    return true;
  }
}

void DirectMPISendNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    runtime.enqueueHost([&] {
      std::lock_guard guard(requestMutex);
      MPI_Startall(requests.size(), requests.data());
    });
    ++progress;
  }
}

void DirectMPISendNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {
  runtime.enqueueHost([&] { MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE); });
}

DirectMPISendNeighborCluster::DirectMPISendNeighborCluster(
    const std::vector<RemoteCluster>& remote,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
    : SendNeighborCluster(cpuExecutor) {
  requests.resize(remote.size());
  status.resize(remote.size());
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

DirectMPISendNeighborCluster::~DirectMPISendNeighborCluster() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
}

bool DirectMPIRecvNeighborCluster::poll() {
  if (requests.size() > 0) {
    std::lock_guard guard(requestMutex);
    int result;
    MPI_Testsome(requests.size(), requests.data(), &result, status.data(), MPI_STATUSES_IGNORE);
    return result == requests.size();
  } else {
    return true;
  }
}

void DirectMPIRecvNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    runtime.enqueueHost([&] {
      std::lock_guard guard(requestMutex);
      MPI_Startall(requests.size(), requests.data());
    });
    ++progress;
  }
}

void DirectMPIRecvNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {
  runtime.enqueueHost([&] { MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE); });
}

DirectMPIRecvNeighborCluster::DirectMPIRecvNeighborCluster(
    const std::vector<RemoteCluster>& remote,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor)
    : RecvNeighborCluster(cpuExecutor) {
  requests.resize(remote.size());
  status.resize(remote.size());
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

DirectMPIRecvNeighborCluster::~DirectMPIRecvNeighborCluster() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
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

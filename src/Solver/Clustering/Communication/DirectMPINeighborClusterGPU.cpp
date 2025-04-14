// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DirectMPINeighborClusterGPU.h"
#include <Memory/MemoryAllocator.h>
#include <Parallel/MPI.h>
#include <Solver/Clustering/Communication/NeighborClusterGPU.h>
#include <mpi.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::solver::clustering::communication {

bool DirectMPISendNeighborClusterGPU::poll() {
  if (requests.size() > 0) {
    std::lock_guard guard(requestMutex);
    int result = 0;
    MPI_Testall(requests.size(), requests.data(), &result, MPI_STATUSES_IGNORE);
    const auto done = result != 0;
    if (done && *progressEnd < progressRestart) {
#pragma omp atomic write
      *progressEnd = progressRestart;
    }
    return done;
  } else {
    return true;
  }
}

void DirectMPISendNeighborClusterGPU::start(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    runtime.enqueueHost([&] {
      std::lock_guard guard(requestMutex);
      ++progressRestart;
      MPI_Startall(requests.size(), requests.data());
    });
    ++progressStart;
  }
}

void DirectMPISendNeighborClusterGPU::stop(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    device::DeviceInstance::getInstance().api->streamWaitMemory(
        runtime.stream(), progressEnd, progressStart);
  }
}

DirectMPISendNeighborClusterGPU::DirectMPISendNeighborClusterGPU(
    const std::vector<RemoteCluster>& remote,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority)
    : SendNeighborCluster(cpuExecutor, priority) {
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
  progressEnd =
      reinterpret_cast<uint32_t*>(device::DeviceInstance::getInstance().api->allocPinnedMem(
          sizeof(uint32_t), device::Destination::CurrentDevice));
  *progressEnd = 0;
}

DirectMPISendNeighborClusterGPU::~DirectMPISendNeighborClusterGPU() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
  device::DeviceInstance::getInstance().api->freePinnedMem(progressEnd);
}

bool DirectMPIRecvNeighborClusterGPU::poll() {
  if (requests.size() > 0) {
    std::lock_guard guard(requestMutex);
    int result = 0;
    MPI_Testall(requests.size(), requests.data(), &result, MPI_STATUSES_IGNORE);
    const auto done = result != 0;
    if (done && *progressEnd < progressRestart) {
#pragma omp atomic write
      *progressEnd = progressRestart;
    }
    return done;
  } else {
    return true;
  }
}

void DirectMPIRecvNeighborClusterGPU::start(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    runtime.enqueueHost([&] {
      std::lock_guard guard(requestMutex);
      ++progressRestart;
      MPI_Startall(requests.size(), requests.data());
    });
    ++progressStart;
  }
}

void DirectMPIRecvNeighborClusterGPU::stop(parallel::runtime::StreamRuntime& runtime) {
  if (requests.size() > 0) {
    device::DeviceInstance::getInstance().api->streamWaitMemory(
        runtime.stream(), progressEnd, progressStart);
  }
}

DirectMPIRecvNeighborClusterGPU::DirectMPIRecvNeighborClusterGPU(
    const std::vector<RemoteCluster>& remote,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority)
    : RecvNeighborCluster(cpuExecutor, priority) {
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
  progressEnd =
      reinterpret_cast<uint32_t*>(device::DeviceInstance::getInstance().api->allocPinnedMem(
          sizeof(uint32_t), device::Destination::CurrentDevice));
  *progressEnd = 0;
}

DirectMPIRecvNeighborClusterGPU::~DirectMPIRecvNeighborClusterGPU() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
  device::DeviceInstance::getInstance().api->freePinnedMem(progressEnd);
}

/*
bool CopyMPISendNeighborClusterGPU::poll() {
    return true;
}

void CopyMPISendNeighborClusterGPU::start(parallel::runtime::StreamRuntime& runtime) {
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

bool CopyMPIRecvNeighborClusterGPU::poll() {
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

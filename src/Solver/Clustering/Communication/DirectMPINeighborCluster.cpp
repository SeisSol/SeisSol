#include "DirectMPINeighborCluster.h"
#include <mpi.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::solver::clustering::communication {

bool DirectMPISendNeighborCluster::poll() {
  std::lock_guard guard(requestMutex);
  int result;
  MPI_Testsome(requests.size(), requests.data(), &result, status.data(), MPI_STATUSES_IGNORE);
  return result == requests.size();
}

void DirectMPISendNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
  runtime.enqueueHost([&] {
    std::lock_guard guard(requestMutex);
    MPI_Startall(requests.size(), requests.data());
  });
}

void DirectMPISendNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {
  runtime.enqueueHost([&] { MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE); });
}

DirectMPISendNeighborCluster::DirectMPISendNeighborCluster(
    const std::vector<RemoteCluster>& remote) {
  requests.resize(remote.size());
  status.resize(remote.size());
  for (std::size_t i = 0; i < remote.size(); ++i) {
    MPI_Send_init(remote[i].data,
                  remote[i].size,
                  remote[i].datatype,
                  remote[i].rank,
                  remote[i].tag,
                  seissol::mpi::MPI::comm(),
                  &requests[i]);
  }
}

DirectMPISendNeighborCluster::~DirectMPISendNeighborCluster() {
  for (auto& request : requests) {
    MPI_Request_free(&request);
  }
}

bool DirectMPIRecvNeighborCluster::poll() {
  std::lock_guard guard(requestMutex);
  int result;
  MPI_Testsome(requests.size(), requests.data(), &result, status.data(), MPI_STATUSES_IGNORE);
  return result == requests.size();
}

void DirectMPIRecvNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
  runtime.enqueueHost([&] {
    std::lock_guard guard(requestMutex);
    MPI_Startall(requests.size(), requests.data());
  });
}

void DirectMPIRecvNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {
  runtime.enqueueHost([&] { MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE); });
}

DirectMPIRecvNeighborCluster::DirectMPIRecvNeighborCluster(
    const std::vector<RemoteCluster>& remote) {
  requests.resize(remote.size());
  status.resize(remote.size());
  for (std::size_t i = 0; i < remote.size(); ++i) {
    MPI_Send_init(remote[i].data,
                  remote[i].size,
                  remote[i].datatype,
                  remote[i].rank,
                  remote[i].tag,
                  seissol::mpi::MPI::comm(),
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

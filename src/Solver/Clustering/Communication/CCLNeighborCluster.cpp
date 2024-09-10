#include "CCLNeighborCluster.h"

#ifdef USE_CCL
#ifdef CCL_NCCL
#include <nccl.h>
using StreamT = cudaStream_t;
#endif
#ifdef CCL_RCCL
#include <rccl/rccl.h>
using StreamT = hipStream_t;
#endif
#include <device.h>

#if REAL_SIZE == 8
constexpr ncclDataType_t CclReal = ncclDataType_t::ncclFloat64;
#elif REAL_SIZE == 4
constexpr ncclDataType_t CclReal = ncclDataType_t::ncclFloat32;
#endif
#endif

namespace seissol::solver::clustering::communication {

bool CCLSendNeighborCluster::poll() { return true; }
void CCLSendNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
#ifdef USE_CCL
  ncclGroupStart();
  for (std::size_t i = 0; i < remote.size(); ++i) {
    const auto& cluster = remote.at(i);
    ncclSend(cluster.data,
             cluster.size,
             CclReal, // TODO(David): use cluster.datatype here instead
             cluster.rank,
             static_cast<ncclComm_t>(comm),
             static_cast<StreamT>(runtime.stream()));
  }
  ncclGroupEnd();
#endif
}
void CCLSendNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {}

CCLSendNeighborCluster::CCLSendNeighborCluster(const std::vector<RemoteCluster>& remote, void* comm)
    : remote(remote), comm(comm) {
#if defined(USE_CCL) && defined(CCL_NCCL)
  for (const auto& cluster : remote) {
    void* handle;
    ncclCommRegister(reinterpret_cast<ncclComm_t>(comm),
                     cluster.data,
                     cluster.size * sizeof(real), // TODO(David): use cluster.datatype here instead
                     &handle);
    memoryHandles.push_back(handle);
  }
#endif
}
CCLSendNeighborCluster::~CCLSendNeighborCluster() {
#if defined(USE_CCL) && defined(CCL_NCCL)
  for (void* handle : memoryHandles) {
    ncclCommDeregister(reinterpret_cast<ncclComm_t>(comm), &handle);
  }
#endif
}

bool CCLRecvNeighborCluster::poll() { return true; }
void CCLRecvNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
#ifdef USE_CCL
  ncclGroupStart();
  for (std::size_t i = 0; i < remote.size(); ++i) {
    const auto& cluster = remote.at(i);
    ncclRecv(cluster.data,
             cluster.size,
             CclReal, // TODO(David): use cluster.datatype here instead
             cluster.rank,
             static_cast<ncclComm_t>(comm),
             static_cast<StreamT>(runtime.stream()));
  }
  ncclGroupEnd();
#endif
}
void CCLRecvNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {}

CCLRecvNeighborCluster::CCLRecvNeighborCluster(const std::vector<RemoteCluster>& remote, void* comm)
    : remote(remote), comm(comm) {
#if defined(USE_CCL) && defined(CCL_NCCL)
  for (const auto& cluster : remote) {
    void* handle;
    ncclCommRegister(reinterpret_cast<ncclComm_t>(comm),
                     cluster.data,
                     cluster.size * sizeof(real), // TODO(David): use cluster.datatype here instead
                     &handle);
    memoryHandles.push_back(handle);
  }
#endif
}
CCLRecvNeighborCluster::~CCLRecvNeighborCluster() {
#if defined(USE_CCL) && defined(CCL_NCCL)
  for (void* handle : memoryHandles) {
    ncclCommDeregister(reinterpret_cast<ncclComm_t>(comm), &handle);
  }
#endif
}

} // namespace seissol::solver::clustering::communication

// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CCLNeighborCluster2.h"
#include <Solver/Clustering/Communication/NeighborCluster.h>

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

namespace {
template <typename T>
constexpr ncclDataType_t cclDatatype() {
  if constexpr (std::is_same_v<T, float>) {
    return ncclDataType_t::ncclFloat32;
  } else if constexpr (std::is_same_v<T, double>) {
    return ncclDataType_t::ncclFloat64;
  } else {
    static_assert(sizeof(T) == 0, "CCL datatype not covered.");
    return ncclDataType_t::ncclChar;
  }
}
} // namespace

#if defined(CCL_NCCL) || (defined(CCL_RCCL) && NCCL_VERSION_CODE >= 22005)
#define USE_CCL_REGISTER
#endif

#endif

namespace seissol::solver::clustering::communication {

bool CCLNeighborCluster::poll() { return true; }
void CCLNeighborCluster::start(parallel::runtime::StreamRuntime& runtime) {
  if (remote.size() > 0) {
#ifdef USE_CCL
    runtime.runGraphGeneric(handle, [&](parallel::runtime::StreamRuntime& runtime) {
      ncclGroupStart();
      for (std::size_t i = 0; i < remote.size(); ++i) {
        const auto& cluster = remote[i];
        if (isSend[i]) {
          ncclSend(cluster.data,
                   cluster.size,
                   cclDatatype<real>(), // TODO(David): use cluster.datatype here instead
                   cluster.rank,
                   static_cast<ncclComm_t>(comm),
                   static_cast<StreamT>(runtime.stream()));
        } else {
          ncclRecv(cluster.data,
                   cluster.size,
                   cclDatatype<real>(), // TODO(David): use cluster.datatype here instead
                   cluster.rank,
                   static_cast<ncclComm_t>(comm),
                   static_cast<StreamT>(runtime.stream()));
        }
      }
      ncclGroupEnd();
    });
#endif
  }
}
void CCLNeighborCluster::stop(parallel::runtime::StreamRuntime& runtime) {}

CCLNeighborCluster::CCLNeighborCluster(
    const std::vector<RemoteCluster>& remoteSend,
    const std::vector<RemoteCluster>& remoteRecv,
    void* comm,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority)
    : comm(comm), SendNeighborCluster(cpuExecutor, priority) {
  remote.resize(remoteRecv.size() + remoteSend.size());
  isSend.resize(remote.size());
  std::size_t recvC = 0;
  std::size_t sendC = 0;
  for (std::size_t i = 0; i < remote.size(); ++i) {
    const bool useSend =
        remoteRecv.size() == recvC || remoteRecv[recvC].tag > remoteSend[sendC].tag;
    isSend[i] = useSend;
    if (useSend) {
      remote[i] = remoteSend[sendC];
      ++sendC;
    } else {
      remote[i] = remoteRecv[recvC];
      ++recvC;
    }
  }

#ifdef USE_CCL_REGISTER
  for (const auto& cluster : remote) {
    void* handle = nullptr;
    ncclCommRegister(reinterpret_cast<ncclComm_t>(comm),
                     cluster.data,
                     cluster.size * sizeof(real), // TODO(David): use cluster.datatype here instead
                     &handle);
    memoryHandles.push_back(handle);
  }
#endif
}

CCLNeighborCluster::~CCLNeighborCluster() {
#ifdef USE_CCL_REGISTER
  for (void* handle : memoryHandles) {
    ncclCommDeregister(reinterpret_cast<ncclComm_t>(comm), &handle);
  }
#endif
}

} // namespace seissol::solver::clustering::communication

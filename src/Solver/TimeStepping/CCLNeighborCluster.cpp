// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "CCLNeighborCluster.h"
#include <Parallel/MPI.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/TimeStepping/AbstractTimeCluster.h>
#include <Solver/TimeStepping/HaloCommunication.h>
#include <memory>

#ifdef USE_CCL
#ifdef CCL_NCCL
#include <nccl.h>
using StreamT = cudaStream_t;
#endif
#ifdef CCL_RCCL
#include <rccl/rccl.h>
using StreamT = hipStream_t;
#endif

#include "utils/logger.h"
#include <device.h>

namespace {
constexpr ncclDataType_t cclDatatype(RealType type) {
  switch (type) {
  case RealType::F32:
    return ncclDataType_t::ncclFloat32;
  case RealType::F64:
    return ncclDataType_t::ncclFloat64;
  default:
    logError() << "Unsupported CCL datatype";
    return ncclDataType_t::ncclChar;
  }
}
} // namespace

#if defined(CCL_NCCL) || (defined(CCL_RCCL) && NCCL_VERSION_CODE >= 22005)
#define USE_CCL_REGISTER
#endif

#endif

namespace {
std::vector<void*> createComms(std::size_t count) {
#ifdef USE_CCL
  // cf. partially https://docs.nvidia.com/deeplearning/nccl/user-guide/docs/examples.html
  std::vector<ncclUniqueId> cclIds(count);
  std::vector<void*> comms(count);
  if (seissol::MPI::mpi.rank() == 0) {
    for (std::size_t i = 0; i < count; ++i) {
      ncclGetUniqueId(&cclIds[i]);
    }
  }
  MPI_Bcast(
      cclIds.data(), sizeof(ncclUniqueId) * cclIds.size(), MPI_BYTE, 0, seissol::MPI::mpi.comm());
  MPI_Barrier(seissol::MPI::mpi.comm());
  for (std::size_t i = 0; i < count; ++i) {
    ncclComm_t preComm = NCCL_COMM_NULL;
    ncclCommInitRank(&preComm, seissol::MPI::mpi.size(), cclIds[i], seissol::MPI::mpi.rank());
    comms[i] = static_cast<void*>(preComm);
  }
  return comms;
#else
  return {};
#endif
}
} // namespace

namespace seissol::time_stepping {

bool CCLNeighborCluster::mayCorrect() { return stream.test() && AbstractTimeCluster::mayCorrect(); }

void CCLNeighborCluster::handleAdvancedPredictionTimeMessage(
    const NeighborCluster& /*neighborCluster*/) {
  assert(stream.test());
  if (remote.size() > 0) {
#ifdef USE_CCL
    stream->runGraphGeneric(handle, [&](parallel::runtime::StreamRuntime& runtime) {
      ncclGroupStart();
      for (std::size_t i = 0; i < remote.size(); ++i) {
        const auto& cluster = remote[i];
        if (isSend[i]) {
          ncclSend(cluster.data,
                   cluster.size,
                   cclDatatype(cluster.datatype),
                   cluster.rank,
                   static_cast<ncclComm_t>(comm),
                   static_cast<StreamT>(runtime.stream()));
        } else {
          ncclRecv(cluster.data,
                   cluster.size,
                   cclDatatype(cluster.datatype),
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

void CCLNeighborCluster::handleAdvancedCorrectionTimeMessage(
    const NeighborCluster& neighborCluster) {}

void CCLNeighborCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logError() << "CCL timeout. More info TBD.";
}

void CCLNeighborCluster::start() {}
void CCLNeighborCluster::predict() {}
void CCLNeighborCluster::correct() {}

CCLNeighborCluster::CCLNeighborCluster(double maxTimeStepSize,
                                       int timeStepRate,
                                       int globalTimeClusterId,
                                       int otherGlobalTimeClusterId,
                                       const seissol::solver::HaloCommunication& meshStructure,
                                       bool persistent)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host) {
  comm = createComms(1)[0];

  // TODO: make session-specific
  static std::shared_ptr<parallel::runtime::StreamRuntime> rt{nullptr};
  if (!rt) {
    rt = std::make_shared<parallel::runtime::StreamRuntime>(0);
  }
  stream = rt;

  const auto& local = meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId);

  remote.resize(local.copy.size() + local.ghost.size());
  isSend.resize(remote.size());

  std::vector<std::pair<solver::RemoteCluster, bool>> temp(remote.size());
  for (std::size_t i = 0; i < local.copy.size(); ++i) {
    temp[i].first = local.copy[i];
    temp[i].second = true;
  }
  for (std::size_t i = 0; i < local.ghost.size(); ++i) {
    temp[i + local.copy.size()].first = local.ghost[i];
    temp[i + local.copy.size()].second = false;
  }

  const auto ownRank = seissol::MPI::mpi.rank();

  std::sort(temp.begin(), temp.end(), [&](const auto& a, const auto& b) {
    const auto& rA = a.first;
    const auto& rB = b.first;

    if (rA.tag == rB.tag) {
      // "hypercube" ordering; cf. e.g.
      // https://github.com/NVIDIA/nccl-tests/blob/abc46770a98777a9fd1b072adcf8becb76bfe125/src/hypercube.cu#L60-L67

      const auto rmA = rA.rank ^ ownRank;
      const auto rmB = rB.rank ^ ownRank;
      return rmA < rmB;
    }

    return rA.tag < rB.tag;
  });

  for (std::size_t i = 0; i < temp.size(); ++i) {
    remote[i] = temp[i].first;
    isSend[i] = temp[i].second;
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

} // namespace seissol::time_stepping

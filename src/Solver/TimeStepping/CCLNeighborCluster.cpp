// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "CCLNeighborCluster.h"

#include "Parallel/MPI.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <memory>

#ifdef USE_CCL

#ifdef CCL_NCCL
#include <nccl.h>
using StreamT = cudaStream_t;
#define CCL(name) nccl##name
#define CCLM(name) NCCL_##name
#endif
#ifdef CCL_RCCL
#include <rccl/rccl.h>
using StreamT = hipStream_t;
#define CCL(name) nccl##name
#define CCLM(name) NCCL_##name
#endif
#ifdef CCL_ONECCL
#include <oneapi/ccl.h>
using StreamT = void*;
#define CCL(name) oneccl##name
#define CCLM(name) ONECCL_##name
#endif

#include <Device/device.h>
#include <utils/logger.h>

namespace seissol {
namespace {
constexpr CCL(DataType_t) cclDatatype(RealType type) {
  switch (type) {
  case RealType::F32:
    return CCL(DataType_t)::CCL(Float32);
  case RealType::F64:
    return CCL(DataType_t)::CCL(Float64);
  default:
    logError() << "Unsupported CCL datatype";
    return CCL(DataType_t)::CCL(Char);
  }
}
} // namespace
} // namespace seissol

#if defined(CCL_NCCL) || (defined(CCL_RCCL) && NCCL_VERSION_CODE >= 22005)
#define USE_CCL_REGISTER
#endif

#endif

namespace {
std::vector<void*> createComms([[maybe_unused]] std::size_t count) {
#ifdef USE_CCL
  // cf. partially https://docs.nvidia.com/deeplearning/nccl/user-guide/docs/examples.html
  std::vector<CCL(UniqueId)> cclIds(count);
  std::vector<void*> comms(count);
  if (seissol::Mpi::mpi.rank() == 0) {
    for (std::size_t i = 0; i < count; ++i) {
      CCL(GetUniqueId)(&cclIds[i]);
    }
  }
  MPI_Bcast(
      cclIds.data(), sizeof(CCL(UniqueId)) * cclIds.size(), MPI_BYTE, 0, seissol::Mpi::mpi.comm());
  MPI_Barrier(seissol::Mpi::mpi.comm());
  for (std::size_t i = 0; i < count; ++i) {
    CCL(Comm_t) preComm = CCLM(COMM_NULL);
    CCL(CommInitRank)(&preComm, seissol::Mpi::mpi.size(), cclIds[i], seissol::Mpi::mpi.rank());
    comms[i] = static_cast<void*>(preComm);
  }
  return comms;
#else
  return {};
#endif
}
} // namespace

namespace seissol::time_stepping {

#ifdef ACL_DEVICE
void CCLNeighborCluster::launch(bool send, bool recv) {
  // we need to split between sends and receives for LTS interfaces
  // (i.e. cluster IDs of copy and ghost don't match)
  // otherwise, we can send+recv everything at once

  auto& handle = [&]() -> device::DeviceGraphHandle& {
    if (send && recv) {
      return handle1;
    }
    if (send && !recv) {
      return handle2;
    }
    if (!send && recv) {
      return handle3;
    }
    throw;
  }();
  stream->runGraphGeneric(handle, [&](parallel::runtime::StreamRuntime& runtime) {
#ifdef USE_CCL
    CCL(GroupStart)();
    for (std::size_t i = 0; i < remote.size(); ++i) {
      const auto& cluster = remote[i];
      if (isSend[i]) {
        if (send) {
          CCL(Send)(cluster.data,
                    cluster.size,
                    cclDatatype(cluster.datatype),
                    cluster.rank,
                    static_cast<CCL(Comm_t)>(comm),
                    static_cast<StreamT>(runtime.stream()));
        }
      } else {
        if (recv) {
          CCL(Recv)(cluster.data,
                    cluster.size,
                    cclDatatype(cluster.datatype),
                    cluster.rank,
                    static_cast<CCL(Comm_t)>(comm),
                    static_cast<StreamT>(runtime.stream()));
        }
      }
    }
    CCL(GroupEnd)();
#endif
  });
}
#endif

bool CCLNeighborCluster::mayCorrect() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  return (event == nullptr || device.api->isEventCompleted(event)) &&
         AbstractTimeCluster::mayCorrect();
#endif
  return AbstractTimeCluster::mayCorrect();
}

void CCLNeighborCluster::handleAdvancedPredictionTimeMessage(
    const NeighborCluster& neighborCluster) {
  if (remote.size() > 0) {
    if (globalClusterId == otherGlobalClusterId) {
      launch(true, true);
      event = stream->eventRecord();
    } else {
      launch(true, false);
      if (globalClusterId > otherGlobalClusterId) {
        launch(false, true);
        event = stream->eventRecord();
      }
    }
  }
}

void CCLNeighborCluster::handleAdvancedCorrectionTimeMessage(
    const NeighborCluster& neighborCluster) {

  auto upcomingCorrectionSteps = ct.stepsSinceLastSync;
  if (state == ActorState::Predicted) {
    upcomingCorrectionSteps = ct.nextCorrectionSteps();
  }

  const bool ignoreMessage = upcomingCorrectionSteps >= ct.stepsUntilSync;

  // If we are already at a sync point, we must not post an additional receive, as otherwise start()
  // posts an additional request. This is also true for the last sync point (i.e. end of
  // simulation), as in this case we do not want to have any hanging request.
  if (!ignoreMessage && remote.size() > 0 && globalClusterId < otherGlobalClusterId) {
    launch(false, true);
  }
}

void CCLNeighborCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logError() << "CCL timeout. More info TBD.";
}

void CCLNeighborCluster::start() {
  if (remote.size() > 0) {
    if (globalClusterId < otherGlobalClusterId) {
      launch(false, true);
    }
  }
}
void CCLNeighborCluster::predict() {}
void CCLNeighborCluster::correct() {}

CCLNeighborCluster::CCLNeighborCluster(double maxTimeStepSize,
                                       int timeStepRate,
                                       int globalTimeClusterId,
                                       int otherGlobalTimeClusterId,
                                       const seissol::solver::HaloCommunication& meshStructure,
                                       bool persistent)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      globalClusterId(globalTimeClusterId), otherGlobalClusterId(otherGlobalTimeClusterId) {
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

  const auto ownRank = seissol::Mpi::mpi.rank();

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
    CCL(CommRegister)(reinterpret_cast<CCL(Comm_t)>(comm),
                      cluster.data,
                      cluster.size * sizeOfRealType(cluster.datatype),
                      &handle);
    memoryHandles.push_back(handle);
  }
#endif
}

CCLNeighborCluster::~CCLNeighborCluster() {
#ifdef USE_CCL_REGISTER
  for (void* handle : memoryHandles) {
    CCL(CommDeregister)(reinterpret_cast<CCL(Comm_t)>(comm), &handle);
  }
#endif
}

std::string CCLNeighborCluster::description() const {
  return "comm-" + std::to_string(globalClusterId) + "-" + std::to_string(otherGlobalClusterId);
}

} // namespace seissol::time_stepping

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CclExchangeScheduler.h"

#include "Solver/TimeStepping/ExchangeScheduler.h"

#include <cstddef>
#include <utils/logger.h>

#if defined(ACL_DEVICE) && defined(USE_CCL)

#include "Kernels/Precision.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <Device/device.h>
#include <algorithm>
#include <mpi.h>
#include <utility>
#include <vector>

// the three libraries offer the same functionality under different prefixes
#ifdef CCL_NCCL
#include <nccl.h>
using StreamT = cudaStream_t;
#define CCL(name) nccl##name
#define CCLM(name) NCCL_##name
#if NCCL_VERSION_CODE >= NCCL_VERSION(2, 19, 0)
#define USE_CCL_REGISTER
#endif
#endif

#ifdef CCL_RCCL
#include <rccl/rccl.h>
using StreamT = hipStream_t;
#define CCL(name) nccl##name
#define CCLM(name) NCCL_##name
#if NCCL_VERSION_CODE >= 22005
#define USE_CCL_REGISTER
#endif
#endif

#ifdef CCL_ONECCL
// the C API (https://github.com/uxlfoundation/oneCCL/blob/rfcs/rfcs/20240806-c-api/README.md)
#include <oneapi/ccl.h>
using StreamT = void*;
#define CCL(name) oneccl##name
#define CCLM(name) ONECCL_##name
#endif

namespace seissol::time_stepping {

namespace {

void check(CCL(Result_t) result, const char* call) {
  if (result != CCL(Success)) {
    logError() << "The CCL call" << call << "failed with the error code"
               << static_cast<int>(result);
  }
}

CCL(DataType_t) datatype(RealType type) {
  switch (type) {
  case RealType::F32:
    return CCL(Float32);
  case RealType::F64:
    return CCL(Float64);
  default:
    logError() << "The halo exchange with CCL supports only single and double precision.";
    return CCL(Char);
  }
}

} // namespace

CclExchangeScheduler::CclExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : StreamExchangeScheduler(clusterCount, order),
      communicators_(order == LaunchOrder::Global ? 1 : clusterCount * clusterCount, nullptr) {
  const auto slots = usedSlots();

  std::vector<CCL(UniqueId)> ids(slots.size());
  if (Mpi::mpi.rank() == 0) {
    for (auto& id : ids) {
      check(CCL(GetUniqueId)(&id), "GetUniqueId");
    }
  }
  MPI_Bcast(ids.data(),
            static_cast<int>(sizeof(CCL(UniqueId)) * ids.size()),
            MPI_BYTE,
            0,
            Mpi::mpi.comm());

  for (std::size_t i = 0; i < slots.size(); ++i) {
    CCL(Comm_t) communicator = CCLM(COMM_NULL);
    check(CCL(CommInitRank)(&communicator, Mpi::mpi.size(), ids[i], Mpi::mpi.rank()),
          "CommInitRank");
    communicators_[slots[i]] = static_cast<void*>(communicator);
  }
}

CclExchangeScheduler::~CclExchangeScheduler() {
  synchronize();
#ifdef USE_CCL_REGISTER
  for (const auto& [communicator, handle] : registrations_) {
    check(CCL(CommDeregister)(static_cast<CCL(Comm_t)>(communicator), handle), "CommDeregister");
  }
#endif
  for (auto* communicator : communicators_) {
    if (communicator != nullptr) {
      check(CCL(CommDestroy)(static_cast<CCL(Comm_t)>(communicator)), "CommDestroy");
    }
  }
}

void CclExchangeScheduler::added([[maybe_unused]] const ScheduledTransport& transport) {
#ifdef USE_CCL_REGISTER
  // the copy regions go out in the direction towards the other cluster, the ghost regions come in
  // from it
  const auto registerRegions =
      [&](const std::vector<solver::RemoteCluster>& regions, std::size_t from, std::size_t to) {
        auto* communicator = communicators_[slot(from, to)];
        for (const auto& region : regions) {
          void* handle = nullptr;
          check(CCL(CommRegister)(static_cast<CCL(Comm_t)>(communicator),
                                  region.data,
                                  region.size * sizeOfRealType(region.datatype),
                                  &handle),
                "CommRegister");
          registrations_.emplace_back(communicator, handle);
        }
      };
  registerRegions(transport.regions().copy, transport.cluster(), transport.otherCluster());
  registerRegions(transport.regions().ghost, transport.otherCluster(), transport.cluster());
#endif
}

void CclExchangeScheduler::enqueueGroup(std::size_t slot,
                                        std::size_t /*from*/,
                                        std::size_t /*to*/,
                                        std::size_t /*exchange*/,
                                        const ScheduledTransport* sender,
                                        const ScheduledTransport* receiver) {
  auto* communicator = static_cast<CCL(Comm_t)>(communicators_[slot]);
  auto* nativeStream = static_cast<StreamT>(stream(slot));

  check(CCL(GroupStart)(), "GroupStart");
  if (sender != nullptr) {
    for (const auto& region : sender->regions().copy) {
      check(CCL(Send)(region.data,
                      region.size,
                      datatype(region.datatype),
                      region.rank,
                      communicator,
                      nativeStream),
            "Send");
    }
  }
  if (receiver != nullptr) {
    for (const auto& region : receiver->regions().ghost) {
      check(CCL(Recv)(region.data,
                      region.size,
                      datatype(region.datatype),
                      region.rank,
                      communicator,
                      nativeStream),
            "Recv");
    }
  }
  check(CCL(GroupEnd)(), "GroupEnd");
}

} // namespace seissol::time_stepping

#else

namespace seissol::time_stepping {

CclExchangeScheduler::CclExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : StreamExchangeScheduler(clusterCount, order) {
  logError() << "This build of SeisSol does not support exchanging the halo data with CCL.";
}

CclExchangeScheduler::~CclExchangeScheduler() = default;

void CclExchangeScheduler::added(const ScheduledTransport& /*transport*/) {}

void CclExchangeScheduler::enqueueGroup(std::size_t /*slot*/,
                                        std::size_t /*from*/,
                                        std::size_t /*to*/,
                                        std::size_t /*exchange*/,
                                        const ScheduledTransport* /*sender*/,
                                        const ScheduledTransport* /*receiver*/) {}

} // namespace seissol::time_stepping

#endif

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CclExchangeScheduler.h"

#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/StreamExchangeScheduler.h"

#include <cstddef>
#include <utils/logger.h>

#if defined(ACL_DEVICE) && defined(USE_CCL)

#include "Common/Real.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"

#include <Device/device.h>
#include <algorithm>
#include <mpi.h>
#include <stdlib.h>
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

namespace seissol::solver {

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

CclExchangeScheduler::CclExchangeScheduler(std::size_t clusterCount,
                                           LaunchOrder order,
                                           std::size_t configCount)
    : StreamExchangeScheduler(clusterCount, order, configCount),
      communicators_(order == LaunchOrder::Global ? 1 : clusterCount * clusterCount, nullptr),
      senders_(communicators_.size(), nullptr), receivers_(communicators_.size(), nullptr) {
#if defined(CCL_NCCL) || defined(CCL_RCCL)
  // After a replay of a recorded super-timestep, the next exchanges are enqueued right away, while
  // those of the replay are still outstanding on the same communicator. That needs the support for
  // mixing graph and direct launches, which NCCL has on by default, but RCCL (7.2) off: its proxy
  // then sends the later messages first, and the processes take them for the wrong operations.
  // A setting of the user is kept. (Only one thread runs while the schedulers are set up.)
  // NOLINTNEXTLINE(concurrency-mt-unsafe)
  setenv("NCCL_GRAPH_MIXING_SUPPORT", "1", 0);
#endif
  // NCCL keeps its state of a capture by its id, which CUDA gives all captures into one graph
  // alike; so the groups of a recorded super-timestep get captured into graphs of their own
  isolateGroups();

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

void CclExchangeScheduler::added(const ScheduledTransport& transport) {
  senders_[slot(transport.cluster(), transport.otherCluster())] = &transport;
  receivers_[slot(transport.otherCluster(), transport.cluster())] = &transport;
#ifdef USE_CCL_REGISTER
  // the copy regions go out in the direction towards the other cluster, the ghost regions come in
  // from it
  const auto registerRegions =
      [&](const std::vector<RemoteCluster>& regions, std::size_t from, std::size_t to) {
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

void CclExchangeScheduler::checkAsyncErrors() const {
#if defined(CCL_NCCL) || defined(CCL_RCCL)
  // e.g. a failed connection, which leaves the operations of the device waiting otherwise
  for (auto* communicator : communicators_) {
    if (communicator == nullptr) {
      continue;
    }
    auto asyncError = CCL(Success);
    const auto result = CCL(CommGetAsyncError)(static_cast<CCL(Comm_t)>(communicator), &asyncError);
    if (result == CCL(Success) && asyncError != CCL(Success) && asyncError != CCL(InProgress)) {
      logError() << "The halo exchange with CCL failed asynchronously:"
                 << CCL(GetErrorString)(asyncError) << "("
                 << CCL(GetLastError)(static_cast<CCL(Comm_t)>(communicator)) << ")";
    }
  }
#endif
}

void CclExchangeScheduler::prepare() {
  if (launchOrder() != LaunchOrder::PerDirection) {
    return;
  }
  // The first group towards a peer on a communicator blocks the host until the peer has joined it
  // there, to set up the connection. With a communicator per direction, two processes could each
  // launch the first group of a different direction and wait for each other forever: e.g. if each
  // of them only receives in one of the two directions between two clusters, it launches that
  // group right at the start. So all communicators get their connections here, one after the other
  // in the same order on all processes, with one group of each direction as it comes later. Its
  // data does not matter; the first exchange overwrites the ghost layers before anything reads
  // them.
  for (const auto current : usedSlots()) {
    enqueueGroup(current, 0, 0, 0, senders_[current], receivers_[current]);
  }
  synchronize();
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

} // namespace seissol::solver

#else

namespace seissol::solver {

CclExchangeScheduler::CclExchangeScheduler(std::size_t clusterCount,
                                           LaunchOrder order,
                                           std::size_t configCount)
    : StreamExchangeScheduler(clusterCount, order, configCount) {
  logError() << "This build of SeisSol does not support exchanging the halo data with CCL.";
}

CclExchangeScheduler::~CclExchangeScheduler() = default;

void CclExchangeScheduler::checkAsyncErrors() const {}

void CclExchangeScheduler::prepare() {}

void CclExchangeScheduler::added(const ScheduledTransport& /*transport*/) {}

void CclExchangeScheduler::enqueueGroup(std::size_t /*slot*/,
                                        std::size_t /*from*/,
                                        std::size_t /*to*/,
                                        std::size_t /*exchange*/,
                                        const ScheduledTransport* /*sender*/,
                                        const ScheduledTransport* /*receiver*/) {}

} // namespace seissol::solver

#endif

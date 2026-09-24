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

device::DeviceInstance& device() { return device::DeviceInstance::getInstance(); }

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
    : ExchangeScheduler(clusterCount, order), clusterCount_(clusterCount) {
  const auto slots = order == LaunchOrder::Global ? 1 : clusterCount * clusterCount;
  communicators_.resize(slots, nullptr);
  streams_.resize(slots, nullptr);

  // only neighboring time clusters exchange data; each ordered pair of them is one direction
  std::vector<std::pair<std::size_t, std::size_t>> directions;
  if (order == LaunchOrder::Global) {
    directions.emplace_back(0, 0);
  } else {
    for (std::size_t from = 0; from < clusterCount; ++from) {
      for (std::size_t to = from == 0 ? 0 : from - 1; to < std::min(from + 2, clusterCount); ++to) {
        directions.emplace_back(from, to);
      }
    }
  }

  std::vector<CCL(UniqueId)> ids(directions.size());
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

  for (std::size_t i = 0; i < directions.size(); ++i) {
    const auto [from, to] = directions[i];
    CCL(Comm_t) communicator = CCLM(COMM_NULL);
    check(CCL(CommInitRank)(&communicator, Mpi::mpi.size(), ids[i], Mpi::mpi.rank()),
          "CommInitRank");
    communicators_[index(from, to)] = static_cast<void*>(communicator);
    streams_[index(from, to)] = device().api->createStream();
  }
}

CclExchangeScheduler::~CclExchangeScheduler() {
  for (auto* stream : streams_) {
    if (stream != nullptr) {
      device().api->syncStreamWithHost(stream);
    }
  }
#ifdef USE_CCL_REGISTER
  for (const auto& [communicator, handle] : registrations_) {
    check(CCL(CommDeregister)(static_cast<CCL(Comm_t)>(communicator), handle), "CommDeregister");
  }
#endif
  for (const auto& [ticket, event] : pendingEvents_) {
    device().api->destroyEvent(event);
  }
  for (auto* communicator : communicators_) {
    if (communicator != nullptr) {
      check(CCL(CommDestroy)(static_cast<CCL(Comm_t)>(communicator)), "CommDestroy");
    }
  }
  for (auto* stream : streams_) {
    if (stream != nullptr) {
      device().api->destroyGenericStream(stream);
    }
  }
}

std::size_t CclExchangeScheduler::index(std::size_t from, std::size_t to) const {
  return launchOrder() == LaunchOrder::Global ? 0 : from * clusterCount_ + to;
}

void CclExchangeScheduler::added([[maybe_unused]] const ScheduledTransport& transport) {
#ifdef USE_CCL_REGISTER
  // the copy regions go out in the direction towards the other cluster, the ghost regions come in
  // from it
  const auto registerRegions =
      [&](const std::vector<solver::RemoteCluster>& regions, std::size_t from, std::size_t to) {
        auto* communicator = communicators_[index(from, to)];
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

ExchangeScheduler::Ticket CclExchangeScheduler::launch(std::size_t from,
                                                       std::size_t to,
                                                       const ScheduledTransport* sender,
                                                       const ScheduledTransport* receiver) {
  auto* communicator = static_cast<CCL(Comm_t)>(communicators_[index(from, to)]);
  auto* stream = streams_[index(from, to)];
  if (communicator == nullptr) {
    logError() << "There is no CCL communicator for the halo exchange from cluster" << from
               << "to cluster" << to;
  }

  check(CCL(GroupStart)(), "GroupStart");
  if (sender != nullptr) {
    for (const auto& region : sender->regions().copy) {
      check(CCL(Send)(region.data,
                      region.size,
                      datatype(region.datatype),
                      region.rank,
                      communicator,
                      static_cast<StreamT>(stream)),
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
                      static_cast<StreamT>(stream)),
            "Recv");
    }
  }
  check(CCL(GroupEnd)(), "GroupEnd");

  auto* event = device().api->createEvent();
  device().api->recordEventOnStream(event, stream);
  const auto ticket = nextTicket_++;
  pendingEvents_[ticket] = event;
  return ticket;
}

bool CclExchangeScheduler::completed(Ticket ticket) {
  const auto pending = pendingEvents_.find(ticket);
  if (pending == pendingEvents_.end()) {
    // tickets are handed out in increasing order; only completed ones are forgotten
    return ticket < nextTicket_;
  }
  if (device().api->isEventCompleted(pending->second)) {
    device().api->destroyEvent(pending->second);
    pendingEvents_.erase(pending);
    return true;
  }
  return false;
}

} // namespace seissol::time_stepping

#else

namespace seissol::time_stepping {

CclExchangeScheduler::CclExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : ExchangeScheduler(clusterCount, order), clusterCount_(clusterCount) {
  logError() << "This build of SeisSol does not support exchanging the halo data with CCL.";
}

CclExchangeScheduler::~CclExchangeScheduler() = default;

void CclExchangeScheduler::added(const ScheduledTransport& /*transport*/) {}

ExchangeScheduler::Ticket CclExchangeScheduler::launch(std::size_t /*from*/,
                                                       std::size_t /*to*/,
                                                       const ScheduledTransport* /*sender*/,
                                                       const ScheduledTransport* /*receiver*/) {
  return 0;
}

bool CclExchangeScheduler::completed(Ticket /*ticket*/) { return true; }

std::size_t CclExchangeScheduler::index(std::size_t from, std::size_t to) const {
  return from * clusterCount_ + to;
}

} // namespace seissol::time_stepping

#endif

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ActorMpiExchangeScheduler.h"

#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/StreamExchangeScheduler.h"

#include <cstddef>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Common/Real.h"
#include "Kernels/DeviceAux/ExchangeEngineAux.h"
#include "Kernels/DeviceAux/HostCopyAux.h"
#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/HostActor.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/Mpi/MpiHaloTransport.h"

#include <Device/device.h>
#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <utils/env.h>
#include <vector>

namespace seissol::solver {

namespace {
using parallel::runtime::Event;
using parallel::runtime::HostActor;
using parallel::runtime::HostCounter;
namespace hostcopy = kernels::device::aux::hostcopy;
namespace exchangeengine = kernels::device::aux::exchangeengine;

// the counters of a transport in host memory
enum Counter : std::size_t {
  // the sends whose data MPI may take (by the device): copied into the host buffer, or without
  // staging, there in the copy regions
  CopiedOut,
  // the sends that MPI has completed (by the host)
  Sent,
  // the receives that MPI has completed (by the host)
  Received,
  // the receives whose ghost regions the device is done with (by the device): copied out of the
  // host buffer, or without staging, ready to be overwritten
  CopiedIn,
  CounterCount
};

// where the regions start in the host buffers: for the vector loads of the copies
constexpr std::size_t RegionAlignment = 256;

device::DeviceInstance& deviceInstance() { return device::DeviceInstance::instance(); }

std::size_t bytes(const RemoteCluster& region) {
  return region.size * sizeOfRealType(region.datatype);
}

std::vector<std::size_t> offsets(const std::vector<RemoteCluster>& regions) {
  std::vector<std::size_t> offsets;
  std::size_t total = 0;
  for (const auto& region : regions) {
    offsets.push_back(total);
    total += (bytes(region) + RegionAlignment - 1) / RegionAlignment * RegionAlignment;
  }
  // the size of the buffer last
  offsets.push_back(total);
  return offsets;
}

/// the regions with their data in `buffer` instead, at `offsets`
std::vector<RemoteCluster> inBuffer(std::vector<RemoteCluster> regions,
                                    void* buffer,
                                    const std::vector<std::size_t>& offsets) {
  for (std::size_t i = 0; i < regions.size(); ++i) {
    regions[i].data = static_cast<std::uint8_t*>(buffer) + offsets[i];
  }
  return regions;
}

void advance(std::uint32_t* counter, std::uint32_t value) {
  // the device reads it while it waits; what MPI has written is there before
  __atomic_store_n(counter, value, __ATOMIC_RELEASE);
}
/// pinned host memory of the given size, if `staged`
void* allocate(bool staged, std::size_t bytes) {
  return staged ? deviceInstance().api().allocPinnedMem(std::max<std::size_t>(bytes, 1)) : nullptr;
}
} // namespace

struct ActorMpiExchangeScheduler::Staging {
  Staging(const RemoteClusterPair& regions, bool persistent, bool staged)
      : regions(regions), copyOffsets(offsets(regions.copy)), ghostOffsets(offsets(regions.ghost)),
        hostCopy(allocate(staged, copyOffsets.back())),
        hostGhost(allocate(staged, ghostOffsets.back())),
        mpi(staged ? RemoteClusterPair{inBuffer(regions.copy, hostCopy, copyOffsets),
                                       inBuffer(regions.ghost, hostGhost, ghostOffsets)}
                   : regions,
            persistent) {
    auto& api = deviceInstance().api();
    counters = static_cast<std::uint32_t*>(api.allocPinnedMem(
        CounterCount * sizeof(std::uint32_t), false, device::Destination::CurrentDevice));
    std::fill_n(counters, CounterCount, 0);
    countersOnDevice = static_cast<std::uint32_t*>(api.devicePointer(counters));
    if (staged) {
      hostCopyOnDevice = static_cast<std::uint8_t*>(api.devicePointer(hostCopy));
      hostGhostOnDevice = static_cast<std::uint8_t*>(api.devicePointer(hostGhost));
    }
    // MPI goes on while the thread that enqueues the work of the device may block in it
    sendQueue.setOnProgressThread(true);
    receiveQueue.setOnProgressThread(true);
  }

  ~Staging() {
    sendQueue.drain();
    receiveQueue.drain();
    mpi.finalize();
    auto& api = deviceInstance().api();
    if (hostCopy != nullptr) {
      api.freePinnedMem(hostCopy);
      api.freePinnedMem(hostGhost);
    }
    api.freePinnedMem(counters);
  }

  Staging(const Staging&) = delete;
  Staging(Staging&&) = delete;
  Staging& operator=(const Staging&) = delete;
  Staging& operator=(Staging&&) = delete;

  // on the device
  RemoteClusterPair regions;
  // where each region starts in its host buffer; the size of the buffer last
  std::vector<std::size_t> copyOffsets;
  std::vector<std::size_t> ghostOffsets;
  // the copy regions and the ghost regions, each one after the other in pinned host memory; and
  // their addresses on the device; none without staging
  void* hostCopy{nullptr};
  void* hostGhost{nullptr};
  std::uint8_t* hostCopyOnDevice{nullptr};
  std::uint8_t* hostGhostOnDevice{nullptr};
  // the counters in pinned host memory, and their address on the device
  std::uint32_t* counters{nullptr};
  std::uint32_t* countersOnDevice{nullptr};
  std::uint32_t sends{0};
  std::uint32_t receives{0};
  // the position among the transports of the engine
  std::size_t index{0};
  // the MPI operations on the host buffers, or without staging, on the regions
  MpiHaloTransport mpi;
  HostActor::Queue sendQueue;
  HostActor::Queue receiveQueue;
};

ActorMpiExchangeScheduler::ActorMpiExchangeScheduler(std::size_t clusterCount,
                                                     LaunchOrder order,
                                                     std::size_t configCount,
                                                     std::size_t lanes,
                                                     bool persistent)
    : StreamExchangeScheduler(clusterCount, order, configCount, lanes), persistent_(persistent) {
  utils::Env env("SEISSOL_");
  const auto wait = env.get<std::string>("ACTOR_MPI_WAIT", "stream");
  if (wait != "stream" && wait != "kernel") {
    logWarning() << "SEISSOL_ACTOR_MPI_WAIT takes stream or kernel; the streams wait as the device"
                 << "offers it.";
  }
  waitKernel_ = wait == "kernel";
  if (waitKernel_) {
    logInfo() << "The streams of the halo exchange wait for the host with a kernel that polls"
              << "(SEISSOL_ACTOR_MPI_WAIT=kernel).";
  }
  staged_ = env.get<bool>("ACTOR_MPI_STAGING", true);
  if (!staged_) {
    logInfo() << "MPI sends and receives the halo data in device memory"
              << "(SEISSOL_ACTOR_MPI_STAGING=0); it must not call the device runtime there.";
  }
  engine_ = env.get<bool>("ACTOR_MPI_ENGINE", false);
  engineBlocks_ = env.get<unsigned>("ACTOR_MPI_ENGINE_BLOCKS", 4);
  engineFromEpoch_ = env.get<unsigned>("ACTOR_MPI_ENGINE_FROM_EPOCH", 0);
  if (engine_ && engineFromEpoch_ > 0) {
    logInfo() << "The exchange engine starts with epoch" << engineFromEpoch_ + 1
              << "(SEISSOL_ACTOR_MPI_ENGINE_FROM_EPOCH).";
  }
  if (engine_ && !staged_) {
    logWarning()
        << "The exchange engine (SEISSOL_ACTOR_MPI_ENGINE=1) needs the staging through host"
        << "memory; the exchange goes without it.";
    engine_ = false;
  }
  if (engine_) {
    logInfo() << "A kernel of" << engineBlocks_
              << "blocks takes up the copies of the halo exchange from and to host memory in the"
              << "order they get ready, during the time stepping (SEISSOL_ACTOR_MPI_ENGINE=1).";
  }
  if (order == LaunchOrder::PerDirection) {
    logWarning()
        << "With the actor-mpi transfer mode, the stream of each direction may wait for the"
        << "host; together, these streams may block all hardware queues of the device.";
  }
}

ActorMpiExchangeScheduler::~ActorMpiExchangeScheduler() {
  // (all MPI operations have completed by now)
  synchronize();
  if (engineRuntime_ != nullptr) {
    quiesce();
    engineRuntime_->wait();
    engineRuntime_->dispose();
  }
  auto& api = deviceInstance().api();
  for (auto* table : engineTables_) {
    api.freeGlobMem(table);
  }
  if (engineTransports_ != nullptr) {
    api.freeGlobMem(engineTransports_);
    api.freePinnedMem(engineControl_);
  }
  staging_.clear();
}

void ActorMpiExchangeScheduler::waitFor(void* stream, std::uint32_t* counter, std::uint32_t value) {
  if (waitKernel_) {
    hostcopy::waitHost(
        static_cast<const std::uint32_t*>(deviceInstance().api().devicePointer(counter)),
        value,
        stream);
  } else {
    deviceInstance().api().streamWaitMemory(stream, counter, value);
  }
}

void ActorMpiExchangeScheduler::added(const ScheduledTransport& transport) {
  auto staging = std::make_unique<Staging>(transport.regions(), persistent_, staged_);
  staging->index = staging_.size();
  staging_[&transport] = std::move(staging);
}

void ActorMpiExchangeScheduler::prepare() {
  if (!engine_) {
    return;
  }
  auto& api = deviceInstance().api();
  const auto upload = [&](const std::vector<exchangeengine::Region>& regions) {
    if (regions.empty()) {
      return static_cast<exchangeengine::Region*>(nullptr);
    }
    const auto size = regions.size() * sizeof(exchangeengine::Region);
    auto* table = static_cast<exchangeengine::Region*>(api.allocGlobMem(size));
    api.copyTo(table, regions.data(), size);
    engineTables_.push_back(table);
    return table;
  };
  engineTable_.resize(staging_.size());
  for (const auto& [transport, staging] : staging_) {
    std::vector<exchangeengine::Region> out;
    for (std::size_t region = 0; region < staging->regions.copy.size(); ++region) {
      out.push_back({staging->hostCopyOnDevice + staging->copyOffsets[region],
                     staging->regions.copy[region].data,
                     bytes(staging->regions.copy[region])});
    }
    std::vector<exchangeengine::Region> in;
    for (std::size_t region = 0; region < staging->regions.ghost.size(); ++region) {
      in.push_back({staging->regions.ghost[region].data,
                    staging->hostGhostOnDevice + staging->ghostOffsets[region],
                    bytes(staging->regions.ghost[region])});
    }
    auto& entry = engineTable_[staging->index];
    entry.out = upload(out);
    entry.outCount = static_cast<std::uint32_t>(out.size());
    entry.in = upload(in);
    entry.inCount = static_cast<std::uint32_t>(in.size());
    entry.copiedOut = &staging->countersOnDevice[CopiedOut];
    entry.sent = &staging->countersOnDevice[Sent];
    entry.received = &staging->countersOnDevice[Received];
    entry.copiedIn = &staging->countersOnDevice[CopiedIn];
  }
  engineTransports_ = static_cast<exchangeengine::Transport*>(api.allocGlobMem(
      std::max<std::size_t>(engineTable_.size(), 1) * sizeof(exchangeengine::Transport)));
  const auto controlSize = 1 + 2 * engineTable_.size();
  engineControl_ = static_cast<std::uint32_t*>(api.allocPinnedMem(
      controlSize * sizeof(std::uint32_t), false, device::Destination::CurrentDevice));
  std::fill_n(engineControl_, controlSize, 0);
  engineControlOnDevice_ = static_cast<std::uint32_t*>(api.devicePointer(engineControl_));
  // on a hardware queue of its own, which it occupies while it runs
  engineRuntime_ = std::make_unique<parallel::runtime::StreamRuntime>(0, 1.0);
}

void ActorMpiExchangeScheduler::startEngine() {
  if (engineRunning_) {
    return;
  }
  // the engine before has stopped, as the device has been synchronized since: it goes on from the
  // exchanges done so far (also those before the time stepping, without the engine)
  for (const auto& [transport, staging] : staging_) {
    auto& entry = engineTable_[staging->index];
    entry.sendRequested = staging->sends;
    entry.sendClaimed = staging->sends;
    entry.receiveRequested = staging->receives;
    entry.receiveClaimed = staging->receives;
  }
  if (!engineTable_.empty()) {
    deviceInstance().api().copyTo(engineTransports_,
                                  engineTable_.data(),
                                  engineTable_.size() * sizeof(exchangeengine::Transport));
  }
  __atomic_store_n(&engineControl_[0], 0U, __ATOMIC_RELEASE);
  exchangeengine::run(engineTransports_,
                      static_cast<std::uint32_t>(engineTable_.size()),
                      engineControlOnDevice_,
                      engineBlocks_,
                      engineRuntime_->stream());
  engineRunning_ = true;
}

void ActorMpiExchangeScheduler::quiesce() {
  if (streamOrdered()) {
    ++epochs_;
  }
  if (!engineRunning_) {
    return;
  }
  // the engine stops once it has done all exchanges launched so far
  const auto count = engineTable_.size();
  for (const auto& [transport, staging] : staging_) {
    engineControl_[1 + staging->index] = staging->sends;
    engineControl_[1 + count + staging->index] = staging->receives;
  }
  __atomic_store_n(&engineControl_[0], 1U, __ATOMIC_RELEASE);
  engineRunning_ = false;
}

ActorEvent ActorMpiExchangeScheduler::transportEvent(const ScheduledTransport& transport) const {
  if (!engine_ || !streamOrdered() || epochs_ < engineFromEpoch_) {
    return ExchangeScheduler::transportEvent(transport);
  }
  // with the engine, the work of the copy layer waits for the copies of its own transport only: the
  // copy regions may change once copied out, and the ghost regions are there once copied in
  const auto& staging = *staging_.at(&transport);
  const HostCounter out{&staging.counters[CopiedOut], staging.sends};
  const HostCounter in{&staging.counters[CopiedIn], staging.receives};
  if (staging.sends == 0) {
    return staging.receives == 0 ? ActorEvent() : ActorEvent(in);
  }
  return staging.receives == 0 ? ActorEvent(out) : ActorEvent(out, in);
}

void ActorMpiExchangeScheduler::enqueueGroup(std::size_t slot,
                                             std::size_t /*from*/,
                                             std::size_t /*to*/,
                                             std::size_t /*exchange*/,
                                             const ScheduledTransport* sender,
                                             const ScheduledTransport* receiver) {
  void* deviceStream = stream(slot);
  if (engine_ && streamOrdered() && epochs_ >= engineFromEpoch_) {
    enqueueForEngine(deviceStream, sender, receiver);
    return;
  }
  if (!staged_) {
    enqueueUnstaged(deviceStream, sender, receiver);
    return;
  }
  // the sends first: the other processes may wait for them in this group
  if (sender != nullptr) {
    auto& staging = *staging_.at(sender);
    const auto send = ++staging.sends;
    if (send > 1) {
      // MPI is done with the copies before
      waitFor(deviceStream, &staging.counters[Sent], send - 1);
    }
    for (std::size_t region = 0; region < staging.regions.copy.size(); ++region) {
      hostcopy::copyToHost(staging.hostCopyOnDevice + staging.copyOffsets[region],
                           staging.regions.copy[region].data,
                           bytes(staging.regions.copy[region]),
                           deviceStream);
    }
    hostcopy::signalHost(&staging.countersOnDevice[CopiedOut], send, deviceStream);
    staging.sendQueue.post(Event(HostCounter{&staging.counters[CopiedOut], send}),
                           [&staging]() { staging.mpi.startSend(); });
    staging.sendQueue.post([&staging]() { return staging.mpi.testSend(); },
                           [&staging, send]() { advance(&staging.counters[Sent], send); });
  }
  if (receiver != nullptr) {
    auto& staging = *staging_.at(receiver);
    const auto receive = ++staging.receives;
    // MPI receives once the device has copied the receive before out of the host buffer
    staging.receiveQueue.post(Event(HostCounter{&staging.counters[CopiedIn], receive - 1}),
                              [&staging]() { staging.mpi.startReceive(); });
    staging.receiveQueue.post(
        [&staging]() { return staging.mpi.testReceive(); },
        [&staging, receive]() { advance(&staging.counters[Received], receive); });
    waitFor(deviceStream, &staging.counters[Received], receive);
    for (std::size_t region = 0; region < staging.regions.ghost.size(); ++region) {
      hostcopy::copyFromHost(staging.regions.ghost[region].data,
                             staging.hostGhostOnDevice + staging.ghostOffsets[region],
                             bytes(staging.regions.ghost[region]),
                             deviceStream);
    }
    hostcopy::signalHost(&staging.countersOnDevice[CopiedIn], receive, deviceStream);
  }
}

void ActorMpiExchangeScheduler::enqueueForEngine(void* deviceStream,
                                                 const ScheduledTransport* sender,
                                                 const ScheduledTransport* receiver) {
  startEngine();
  // the stream only requests the copies, once the copy regions are written and the ghost regions
  // may change; the engine does them once MPI is ready for them
  exchangeengine::Transport* sending = nullptr;
  exchangeengine::Transport* receiving = nullptr;
  std::uint32_t send = 0;
  std::uint32_t receive = 0;
  if (sender != nullptr) {
    auto& staging = *staging_.at(sender);
    send = ++staging.sends;
    sending = engineTransports_ + staging.index;
    staging.sendQueue.post(Event(HostCounter{&staging.counters[CopiedOut], send}),
                           [&staging]() { staging.mpi.startSend(); });
    staging.sendQueue.post([&staging]() { return staging.mpi.testSend(); },
                           [&staging, send]() { advance(&staging.counters[Sent], send); });
  }
  if (receiver != nullptr) {
    auto& staging = *staging_.at(receiver);
    receive = ++staging.receives;
    receiving = engineTransports_ + staging.index;
    staging.receiveQueue.post(Event(HostCounter{&staging.counters[CopiedIn], receive - 1}),
                              [&staging]() { staging.mpi.startReceive(); });
    staging.receiveQueue.post(
        [&staging]() { return staging.mpi.testReceive(); },
        [&staging, receive]() { advance(&staging.counters[Received], receive); });
  }
  exchangeengine::request(sending, send, receiving, receive, deviceStream);
}

void ActorMpiExchangeScheduler::enqueueUnstaged(void* deviceStream,
                                                const ScheduledTransport* sender,
                                                const ScheduledTransport* receiver) {
  // MPI may take the copy regions and overwrite the ghost regions from here on; the stream goes on
  // once it is done with both
  std::uint32_t send = 0;
  std::uint32_t receive = 0;
  if (sender != nullptr) {
    auto& staging = *staging_.at(sender);
    send = ++staging.sends;
    hostcopy::signalHost(&staging.countersOnDevice[CopiedOut], send, deviceStream);
    staging.sendQueue.post(Event(HostCounter{&staging.counters[CopiedOut], send}),
                           [&staging]() { staging.mpi.startSend(); });
    staging.sendQueue.post([&staging]() { return staging.mpi.testSend(); },
                           [&staging, send]() { advance(&staging.counters[Sent], send); });
  }
  if (receiver != nullptr) {
    auto& staging = *staging_.at(receiver);
    receive = ++staging.receives;
    hostcopy::signalHost(&staging.countersOnDevice[CopiedIn], receive, deviceStream);
    staging.receiveQueue.post(Event(HostCounter{&staging.counters[CopiedIn], receive}),
                              [&staging]() { staging.mpi.startReceive(); });
    staging.receiveQueue.post(
        [&staging]() { return staging.mpi.testReceive(); },
        [&staging, receive]() { advance(&staging.counters[Received], receive); });
  }
  if (sender != nullptr) {
    waitFor(deviceStream, &staging_.at(sender)->counters[Sent], send);
  }
  if (receiver != nullptr) {
    waitFor(deviceStream, &staging_.at(receiver)->counters[Received], receive);
  }
}

} // namespace seissol::solver

#else

namespace seissol::solver {

struct ActorMpiExchangeScheduler::Staging {};

ActorMpiExchangeScheduler::ActorMpiExchangeScheduler(std::size_t clusterCount,
                                                     LaunchOrder order,
                                                     std::size_t configCount,
                                                     std::size_t lanes,
                                                     bool /*persistent*/)
    : StreamExchangeScheduler(clusterCount, order, configCount, lanes) {
  logError() << "The actor-mpi transfer mode needs a device.";
}

ActorMpiExchangeScheduler::~ActorMpiExchangeScheduler() = default;

void ActorMpiExchangeScheduler::added(const ScheduledTransport& /*transport*/) {}

void ActorMpiExchangeScheduler::prepare() {}

ActorEvent ActorMpiExchangeScheduler::transportEvent(const ScheduledTransport& transport) const {
  return ExchangeScheduler::transportEvent(transport);
}

void ActorMpiExchangeScheduler::quiesce() {}

void ActorMpiExchangeScheduler::enqueueGroup(std::size_t /*slot*/,
                                             std::size_t /*from*/,
                                             std::size_t /*to*/,
                                             std::size_t /*exchange*/,
                                             const ScheduledTransport* /*sender*/,
                                             const ScheduledTransport* /*receiver*/) {}

} // namespace seissol::solver

#endif

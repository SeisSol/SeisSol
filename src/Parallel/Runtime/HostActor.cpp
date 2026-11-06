// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostActor.h"

#include "Parallel/Runtime/Event.h"
#include "Parallel/Runtime/HostBufferPool.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Parallel/Runtime/Stream.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <thread>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Kernels/DeviceAux/HostCopyAux.h"

#include <Device/device.h>
#endif

namespace seissol::parallel::runtime {

namespace {
// where staged data starts in device memory, and in the copies on the host: for vector loads of
// the work, and the copy kernel
constexpr std::size_t StagingAlignment = 256;
// device memory of a recording comes in chunks of at least this many bytes
constexpr std::size_t ChunkBytes = std::size_t{1} << 20U;

std::size_t aligned(std::size_t bytes) {
  return (bytes + StagingAlignment - 1) / StagingAlignment * StagingAlignment;
}
} // namespace

HostActor::Queue::Queue(std::size_t bytes) : Queue(HostActor::instance(), bytes) {}

HostActor::Queue::Queue(HostActor& actor, std::size_t bytes) : actor_(actor), bytes_(bytes) {
  if (bytes_ > 0) {
    buffers_ = std::make_unique<HostBufferPool>(bytes_);
  }
  actor_.queues_.push_back(this);
}

HostActor::Queue::~Queue() {
  // whatever the work reads or writes may go right after this; it has to be done before
  drain();
  auto& queues = actor_.queues_;
  queues.erase(std::remove(queues.begin(), queues.end(), this), queues.end());
}

void HostActor::Queue::post(Event after, std::function<void()> work) {
  if (after.graph() != nullptr) {
    logError() << "Work on the host cannot wait for a node of a graph that is being built.";
  }
  items_.push_back(Item{std::move(after), std::move(work)});
  ++posted_;
}

HostActor::Staging HostActor::Queue::take(StreamRuntime& runtime) {
  if (actor_.recording_ != nullptr) {
    return Staging{actor_.recording_->stage(bytes_, runtime), nullptr};
  }
  const auto buffer = buffers_->take();
  return Staging{buffer.device, buffer.host};
}

void HostActor::Queue::post(StreamRuntime& runtime,
                            const Staging& staging,
                            std::function<void(const void* data)> work) {
  if (actor_.recording_ != nullptr) {
    actor_.recording_->add(*this, staging.device, std::move(work));
    return;
  }
  post(runtime.eventRecord(), [this, staging, work = std::move(work)]() {
    std::invoke(work, staging.host);
    buffers_->give(HostBufferPool::Buffer{staging.host, staging.device});
  });
}

bool HostActor::Queue::progress() {
  bool handedOver = false;
  while (!items_.empty() && actor_.ready(items_.front().after)) {
    // out of the queue first: the work may post more
    auto task = [this, work = std::move(items_.front().work)]() {
      std::invoke(work);
      done_.fetch_add(1, std::memory_order_release);
    };
    items_.pop_front();
    if (actor_.executor_) {
      actor_.executor_(std::move(task));
    } else {
      // the work changes only what belongs to the queue
      runtime_.action(HostAccess{}, std::move(task));
    }
    handedOver = true;
  }
  return handedOver;
}

void HostActor::Queue::drain() {
  if (!items_.empty() && recordingOuterGraph()) {
    logError() << "Work on the host cannot be waited for while a graph is being recorded.";
  }
  while (!items_.empty()) {
    if (!progress()) {
      std::this_thread::yield();
    }
  }
  // the work handed over to tasks or to the executor
  while (!empty()) {
    HostRuntime::yield();
    std::this_thread::yield();
  }
}

HostActor::HostActor()
    : HostActor([]([[maybe_unused]] const Event& event) {
#ifdef ACL_DEVICE
        return device::DeviceInstance::instance().api().isEventCompleted(event.get());
#else
        return true;
#endif
      }) {
}

HostActor::HostActor(Completed completed) : completed_(std::move(completed)) {}

HostActor& HostActor::instance() {
  static HostActor actor;
  return actor;
}

void HostActor::setExecutor(Executor executor) {
  if (pending() > 0) {
    logError() << "The host actor cannot change where its work runs while some is pending.";
  }
  executor_ = std::move(executor);
}

bool HostActor::progress() {
  // (no event may be queried while a graph is being recorded)
  if (recordingOuterGraph()) {
    return false;
  }
  bool handedOver = false;
  // by index: the work may set up queues
  for (std::size_t i = 0; i < queues_.size(); ++i) {
    handedOver = queues_[i]->progress() || handedOver;
  }
  return handedOver;
}

void HostActor::drain() {
  for (std::size_t i = 0; i < queues_.size(); ++i) {
    queues_[i]->drain();
  }
}

void HostActor::beginRecording(Recording& recording) {
  if (recording_ != nullptr) {
    logError() << "The host actor records into one recording at a time.";
  }
  recording_ = &recording;
}

void HostActor::endRecording() { recording_ = nullptr; }

HostActor::Recording::~Recording() {
#ifdef ACL_DEVICE
  for (const auto& chunk : chunks_) {
    device::DeviceInstance::instance().api().freeGlobMem(chunk.device);
  }
#endif
}

void* HostActor::Recording::stage(std::size_t bytes, StreamRuntime& runtime) {
  const auto size = aligned(std::max<std::size_t>(bytes, 1));
  if (chunks_.empty() || chunks_.back().used + size > chunks_.back().capacity) {
#ifdef ACL_DEVICE
    // no memory may be allocated while the stream of the open node gets captured
    runtime.pauseCapture();
    Chunk chunk;
    chunk.capacity = std::max(size, ChunkBytes);
    chunk.device = device::DeviceInstance::instance().api().allocGlobMem(chunk.capacity);
    chunks_.push_back(chunk);
#else
    logError() << "Only a device records graphs, and data for work on the host in them.";
#endif
  }
  auto& chunk = chunks_.back();
  void* staged = static_cast<std::uint8_t*>(chunk.device) + chunk.used;
  chunk.used += size;
  return staged;
}

void HostActor::Recording::add(Queue& queue,
                               const void* staged,
                               std::function<void(const void* data)> work) {
  Item item;
  item.queue = &queue;
  item.work = std::move(work);
  for (std::size_t chunk = 0; chunk < chunks_.size(); ++chunk) {
    const auto* begin = static_cast<const std::uint8_t*>(chunks_[chunk].device);
    const auto* at = static_cast<const std::uint8_t*>(staged);
    if (at >= begin && at < begin + chunks_[chunk].used) {
      item.chunk = chunk;
      item.offset = static_cast<std::size_t>(at - begin);
    }
  }
  items_.push_back(std::move(item));
}

void HostActor::Recording::launched([[maybe_unused]] StreamRuntime& runtime) {
  if (items_.empty()) {
    return;
  }
#ifdef ACL_DEVICE
  if (copies_ == nullptr) {
    // the chunks one after another, as far as they are used
    std::size_t total = 0;
    for (const auto& chunk : chunks_) {
      copyOffsets_.push_back(total);
      total += aligned(chunk.used);
    }
    copies_ = std::make_shared<HostBufferPool>(total);
  }
  const auto buffer = copies_->take();
  for (std::size_t chunk = 0; chunk < chunks_.size(); ++chunk) {
    // by a kernel, as the copy engines may be busy with a large copy of another stream
    kernels::device::aux::hostcopy::copyToHost(static_cast<std::uint8_t*>(buffer.device) +
                                                   copyOffsets_[chunk],
                                               chunks_[chunk].device,
                                               chunks_[chunk].used,
                                               runtime.stream());
  }
  const auto copied = runtime.eventRecord();
  // the copy goes back once all work that reads it is done (also if that is after the recording)
  const auto release = std::shared_ptr<void>(
      nullptr, [copies = copies_, buffer](void* /*unused*/) { copies->give(buffer); });
  for (const auto& item : items_) {
    const auto* data =
        static_cast<const std::uint8_t*>(buffer.host) + copyOffsets_[item.chunk] + item.offset;
    item.queue->post(copied, [work = item.work, data, release]() { std::invoke(work, data); });
  }
#endif
}

std::size_t HostActor::Recording::bytes() const {
  std::size_t bytes = 0;
  for (const auto& chunk : chunks_) {
    bytes += chunk.capacity;
  }
  return bytes;
}

std::size_t HostActor::pending() const {
  std::size_t pending = 0;
  for (const auto* queue : queues_) {
    pending += queue->pending();
  }
  return pending;
}

bool HostActor::ready(const Event& event) const { return !event || completed_(event); }

} // namespace seissol::parallel::runtime

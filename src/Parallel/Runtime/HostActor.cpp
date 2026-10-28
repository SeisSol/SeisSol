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
#include <mutex>
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
  const std::scoped_lock lock(actor_.queuesMutex_);
  actor_.queues_.push_back(this);
}

HostActor::Queue::~Queue() {
  // whatever the work reads or writes may go right after this; it has to be done before
  drain();
  const std::scoped_lock lock(actor_.queuesMutex_);
  auto& queues = actor_.queues_;
  queues.erase(std::remove(queues.begin(), queues.end(), this), queues.end());
}

void HostActor::Queue::post(Event after, std::function<void()> work) {
  if (after.graph() != nullptr) {
    logError() << "Work on the host can only wait for device events, and for counters.";
  }
  if (onProgressThread_ && after.get() != nullptr) {
    logError() << "The progress thread of the host actor does not ask the device about events.";
  }
  const std::scoped_lock lock(mutex_);
  items_.push_back(Item{std::move(after), {}, std::move(work)});
  posted_.fetch_add(1, std::memory_order_release);
}

void HostActor::Queue::post(std::function<bool()> ready, std::function<void()> work) {
  const std::scoped_lock lock(mutex_);
  items_.push_back(Item{Event(), std::move(ready), std::move(work)});
  posted_.fetch_add(1, std::memory_order_release);
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
  // the progress thread takes its queues while it runs
  if (onProgressThread_ && actor_.progressThreadRuns()) {
    return false;
  }
  return progressHere();
}

bool HostActor::Queue::progressHere() {
  bool handedOver = false;
  while (true) {
    std::function<void()> work;
    {
      // (only one thread takes the work of a queue; the lock is for those that post)
      std::unique_lock lock(mutex_);
      if (items_.empty()) {
        break;
      }
      auto& front = items_.front();
      if (!(front.ready ? std::invoke(front.ready) : actor_.ready(front.after))) {
        break;
      }
      // out of the queue first: the work may post more
      work = std::move(front.work);
      items_.pop_front();
    }
    auto task = [this, work = std::move(work)]() {
      std::invoke(work);
      done_.fetch_add(1, std::memory_order_release);
    };
    if (onProgressThread_) {
      std::invoke(task);
    } else if (actor_.executor_) {
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
  if (!empty() && recordingOuterGraph()) {
    logError() << "Work on the host cannot be waited for while a graph is being recorded.";
  }
  // the work of all queues goes on meanwhile: what the device does before this work may wait for
  // what other queues do (e.g. for the MPI operations of the halo exchange); and that handed over
  // to tasks, to the executor, or to the progress thread
  while (!empty()) {
    if (!actor_.progress()) {
      HostRuntime::yield();
      std::this_thread::yield();
    }
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
  // by index: the work may set up queues (on the thread that polls, which changes the list)
  for (std::size_t i = 0; i < queues_.size(); ++i) {
    handedOver = queues_[i]->progress() || handedOver;
  }
  return handedOver;
}

void HostActor::startProgressThread(std::function<void()> onStart) {
  if (progressThread_.joinable()) {
    return;
  }
  stopProgress_.store(false, std::memory_order_relaxed);
  progressThread_ =
      std::thread([this, onStart = std::move(onStart)]() { runProgressThread(onStart); });
}

void HostActor::stopProgressThread() {
  if (progressThread_.joinable()) {
    stopProgress_.store(true, std::memory_order_release);
    progressThread_.join();
  }
}

HostActor::~HostActor() { stopProgressThread(); }

void HostActor::runProgressThread(const std::function<void()>& onStart) {
  if (onStart) {
    std::invoke(onStart);
  }
  // until stopped, and its queues are empty
  while (true) {
    bool handedOver = false;
    bool pending = false;
    {
      const std::scoped_lock lock(queuesMutex_);
      for (auto* queue : queues_) {
        if (queue->onProgressThread_) {
          handedOver = queue->progressHere() || handedOver;
          pending = pending || !queue->empty();
        }
      }
    }
    if (!pending && stopProgress_.load(std::memory_order_acquire)) {
      break;
    }
    if (!handedOver) {
      std::this_thread::yield();
    }
  }
}

void HostActor::drain() {
  for (std::size_t i = 0; i < queues_.size(); ++i) {
    // (each one goes on with all of them)
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

bool HostActor::ready(const Event& event) const {
  if (event.onCounters()) {
    // in host memory: no need to ask the device
    return std::all_of(
        event.counters().begin(), event.counters().end(), [](const HostCounter& counter) {
          return counter.counter == nullptr ||
                 __atomic_load_n(counter.counter, __ATOMIC_ACQUIRE) >= counter.value;
        });
  }
  return !event || completed_(event);
}

} // namespace seissol::parallel::runtime

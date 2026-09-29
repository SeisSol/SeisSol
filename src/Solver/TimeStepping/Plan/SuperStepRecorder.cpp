// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SuperStepRecorder.h"

#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"

#include <vector>

#ifdef ACL_DEVICE
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/Actor/ActorState.h"

#include <Device/device.h>
#include <cstddef>
#include <utils/logger.h>
#endif

namespace seissol::solver {

#ifdef ACL_DEVICE

namespace {
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::instance(); }

// enough to never re-record an event that someone may still have to wait for
constexpr std::size_t EventCount = 16;
} // namespace

SuperStepRecorder::SuperStepRecorder() {
  stream_ = deviceInstance().api().createStream();
  for (std::size_t i = 0; i < EventCount; ++i) {
    events_.push_back(deviceInstance().api().createEvent());
  }
}

SuperStepRecorder::~SuperStepRecorder() { dispose(); }

void SuperStepRecorder::dispose() {
  if (stream_ != nullptr) {
    deviceInstance().api().syncStreamWithHost(stream_);
    // the graphs own their executable instances; they have to go before the device does
    graphs_.clear();
    recording_.reset();
    waitFor_.clear();
    for (auto* event : events_) {
      deviceInstance().api().destroyEvent(event);
    }
    events_.clear();
    deviceInstance().api().destroyGenericStream(stream_);
    stream_ = nullptr;
  }
}

bool SuperStepRecorder::available() { return deviceInstance().api().isCapableOfGraphCapturing(); }

void* SuperStepRecorder::nextEvent() {
  auto* event = events_[eventIndex_];
  eventIndex_ = (eventIndex_ + 1) % events_.size();
  return event;
}

bool SuperStepRecorder::has(const Key& key) const { return graphs_.find(key) != graphs_.end(); }

void* SuperStepRecorder::lastEvent() const { return lastEvent_; }

void SuperStepRecorder::beginRecording(const std::vector<AbstractTimeCluster*>& clusters,
                                       const std::vector<void*>& streams) {
  beginReplay(clusters, streams);

  recording_ = deviceInstance().api().streamBeginCapture({stream_});
  if (!recording_.isInitialized()) {
    logError() << "Could not start recording a super-timestep.";
  }

  // fork all streams; inside the recording, they may only wait for each other
  auto* fork = nextEvent();
  deviceInstance().api().recordEventOnStream(fork, stream_);
  for (auto* cluster : clusters) {
    cluster->joinEvent(fork);
    cluster->publishEvent(ActorEvent(fork));
  }
  for (auto* stream : streams) {
    deviceInstance().api().syncStreamWithEvent(stream, fork);
  }
  parallel::runtime::recordingOuterGraph() = true;
}

void SuperStepRecorder::endRecording(const Key& key,
                                     const std::vector<AbstractTimeCluster*>& clusters,
                                     const std::vector<void*>& streams) {
  parallel::runtime::recordingOuterGraph() = false;

  // join all streams back
  for (auto* cluster : clusters) {
    const auto join = cluster->markEvent();
    if (join) {
      deviceInstance().api().syncStreamWithEvent(stream_, join.get());
    }
  }
  for (auto* stream : streams) {
    auto* join = nextEvent();
    deviceInstance().api().recordEventOnStream(join, stream);
    deviceInstance().api().syncStreamWithEvent(stream_, join);
  }
  deviceInstance().api().streamEndCapture(recording_);
  graphs_[key] = recording_;
}

void SuperStepRecorder::beginReplay(const std::vector<AbstractTimeCluster*>& clusters,
                                    const std::vector<void*>& streams) {
  waitFor_.clear();
  for (auto* cluster : clusters) {
    waitFor_.push_back(cluster->latestEvent());
  }
  for (auto* stream : streams) {
    auto* latest = nextEvent();
    deviceInstance().api().recordEventOnStream(latest, stream);
    waitFor_.emplace_back(latest);
  }
}

void SuperStepRecorder::replay(const Key& key,
                               const std::vector<AbstractTimeCluster*>& clusters,
                               const std::vector<void*>& streams) {
  for (const auto& event : waitFor_) {
    if (event) {
      deviceInstance().api().syncStreamWithEvent(stream_, event.get());
    }
  }
  waitFor_.clear();
  deviceInstance().api().launchGraph(graphs_.at(key), stream_);

  // everything enqueued from now on comes after the replayed work
  auto* done = nextEvent();
  deviceInstance().api().recordEventOnStream(done, stream_);
  lastEvent_ = done;
  for (auto* cluster : clusters) {
    cluster->joinEvent(done);
    cluster->publishEvent(ActorEvent(done));
  }
  for (auto* stream : streams) {
    deviceInstance().api().syncStreamWithEvent(stream, done);
  }
}

#else

SuperStepRecorder::SuperStepRecorder() = default;

SuperStepRecorder::~SuperStepRecorder() = default;

void SuperStepRecorder::dispose() {}

bool SuperStepRecorder::available() { return false; }

bool SuperStepRecorder::has(const Key& /*key*/) const { return false; }

void* SuperStepRecorder::lastEvent() const { return nullptr; }

void SuperStepRecorder::beginRecording(const std::vector<AbstractTimeCluster*>& /*clusters*/,
                                       const std::vector<void*>& /*streams*/) {}

void SuperStepRecorder::endRecording(const Key& /*key*/,
                                     const std::vector<AbstractTimeCluster*>& /*clusters*/,
                                     const std::vector<void*>& /*streams*/) {}

void SuperStepRecorder::beginReplay(const std::vector<AbstractTimeCluster*>& /*clusters*/,
                                    const std::vector<void*>& /*streams*/) {}

void SuperStepRecorder::replay(const Key& /*key*/,
                               const std::vector<AbstractTimeCluster*>& /*clusters*/,
                               const std::vector<void*>& /*streams*/) {}

#endif

} // namespace seissol::solver

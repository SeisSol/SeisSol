// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SuperStepRecorder.h"

#include "Common/Typedefs.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"

#include <cstddef>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::time_stepping {

#ifdef ACL_DEVICE

namespace {
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::getInstance(); }

// enough to never re-record an event that someone may still have to wait for
constexpr std::size_t EventCount = 4;
} // namespace

SuperStepRecorder::SuperStepRecorder() {
  stream_ = deviceInstance().api->createStream();
  for (std::size_t i = 0; i < EventCount; ++i) {
    events_.push_back(deviceInstance().api->createEvent());
  }
}

SuperStepRecorder::~SuperStepRecorder() { dispose(); }

void SuperStepRecorder::dispose() {
  if (stream_ != nullptr) {
    deviceInstance().api->syncStreamWithHost(stream_);
    for (auto* event : events_) {
      deviceInstance().api->destroyEvent(event);
    }
    events_.clear();
    deviceInstance().api->destroyGenericStream(stream_);
    stream_ = nullptr;
  }
}

bool SuperStepRecorder::available() {
  return Backend != DeviceBackend::Hip && deviceInstance().api->isCapableOfGraphCapturing();
}

void* SuperStepRecorder::nextEvent() {
  auto* event = events_[eventIndex_];
  eventIndex_ = (eventIndex_ + 1) % events_.size();
  return event;
}

bool SuperStepRecorder::has(const Key& key) const { return graphs_.find(key) != graphs_.end(); }

void SuperStepRecorder::beginRecording(const std::vector<AbstractTimeCluster*>& clusters) {
  beginReplay(clusters);

  std::vector<void*> streams{stream_};
  recording_ = deviceInstance().api->streamBeginCapture(streams);

  // fork the streams of all clusters; inside the recording, they may only wait for each other
  auto* fork = nextEvent();
  deviceInstance().api->recordEventOnStream(fork, stream_);
  for (auto* cluster : clusters) {
    cluster->joinEvent(fork);
    cluster->publishEvent(fork);
  }
  parallel::runtime::recordingOuterGraph() = true;
}

void SuperStepRecorder::endRecording(const Key& key,
                                     const std::vector<AbstractTimeCluster*>& clusters) {
  parallel::runtime::recordingOuterGraph() = false;

  // join the streams of all clusters back
  for (auto* cluster : clusters) {
    auto* join = cluster->markEvent();
    deviceInstance().api->syncStreamWithEvent(stream_, join);
  }
  deviceInstance().api->streamEndCapture(recording_);
  graphs_[key] = recording_;
}

void SuperStepRecorder::beginReplay(const std::vector<AbstractTimeCluster*>& clusters) {
  waitFor_.clear();
  for (auto* cluster : clusters) {
    waitFor_.push_back(cluster->latestEvent());
  }
}

void SuperStepRecorder::replay(const Key& key, const std::vector<AbstractTimeCluster*>& clusters) {
  for (auto* event : waitFor_) {
    if (event != nullptr) {
      deviceInstance().api->syncStreamWithEvent(stream_, event);
    }
  }
  deviceInstance().api->launchGraph(graphs_.at(key), stream_);

  // everything the clusters enqueue from now on comes after the replayed work
  auto* done = nextEvent();
  deviceInstance().api->recordEventOnStream(done, stream_);
  for (auto* cluster : clusters) {
    cluster->joinEvent(done);
    cluster->publishEvent(done);
  }
}

#else

SuperStepRecorder::SuperStepRecorder() = default;

SuperStepRecorder::~SuperStepRecorder() = default;

void SuperStepRecorder::dispose() {}

bool SuperStepRecorder::available() { return false; }

bool SuperStepRecorder::has(const Key& /*key*/) const { return false; }

void SuperStepRecorder::beginRecording(const std::vector<AbstractTimeCluster*>& /*clusters*/) {}

void SuperStepRecorder::endRecording(const Key& /*key*/,
                                     const std::vector<AbstractTimeCluster*>& /*clusters*/) {}

void SuperStepRecorder::beginReplay(const std::vector<AbstractTimeCluster*>& /*clusters*/) {}

void SuperStepRecorder::replay(const Key& /*key*/,
                               const std::vector<AbstractTimeCluster*>& /*clusters*/) {}

#endif

} // namespace seissol::time_stepping

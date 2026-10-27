// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SuperStepRecorder.h"

#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"

#include <vector>

#ifdef ACL_DEVICE
#include "Parallel/Helper.h"
#include "Parallel/Runtime/Stream.h"

#include <Device/device.h>
#include <algorithm>
#include <iterator>
#include <utility>
#endif

namespace seissol::solver {

#ifdef ACL_DEVICE

namespace {
device::DeviceInstance& deviceInstance() { return device::DeviceInstance::instance(); }
} // namespace

SuperStepRecorder::SuperStepRecorder(bool branches)
    : stream_(runtime_.stream()), branches_(branches) {}

SuperStepRecorder::~SuperStepRecorder() { dispose(); }

void SuperStepRecorder::dispose() {
  if (stream_ != nullptr) {
    deviceInstance().api().syncStreamWithHost(stream_);
    // the graphs own their executable instances; they have to go before the device does
    graphs_.clear();
    used_.clear();
    launchedOnce_.clear();
    waitFor_.clear();
    lastEvent_ = ActorEvent();
    runtime_.dispose();
    stream_ = nullptr;
  }
}

bool SuperStepRecorder::available() { return useGraphNodes(); }

bool SuperStepRecorder::has(const Key& key) const { return graphs_.find(key) != graphs_.end(); }

ActorEvent SuperStepRecorder::lastEvent() const { return lastEvent_; }

void SuperStepRecorder::beginRecording(const std::vector<AbstractTimeCluster*>& clusters,
                                       const Runtimes& runtimes) {
  beginReplay(runtimes);

  runtime_.beginGraph();
  for (auto* runtime : runtimes) {
    runtime->attach(runtime_, branches_);
  }
  for (auto* cluster : clusters) {
    cluster->publishEvent(ActorEvent());
  }
  parallel::runtime::recordingOuterGraph() = true;
}

void SuperStepRecorder::endRecording(const Key& key) {
  parallel::runtime::recordingOuterGraph() = false;
  graphs_[key] = runtime_.endGraph();
  used_.insert(key);
}

void SuperStepRecorder::endRecordingAndLaunch(const std::vector<AbstractTimeCluster*>& clusters,
                                              const Runtimes& runtimes,
                                              double* deviceTime) {
  parallel::runtime::recordingOuterGraph() = false;
  auto recording = runtime_.endGraph();
  // the graphs launched once before that have run by now go (no event may be queried while
  // recording)
  launchedOnce_.erase(std::remove_if(launchedOnce_.begin(),
                                     launchedOnce_.end(),
                                     [](const auto& launched) {
                                       return deviceInstance().api().isEventCompleted(
                                           launched.second.get());
                                     }),
                      launchedOnce_.end());
  launch(recording, clusters, runtimes, deviceTime);
  launchedOnce_.emplace_back(std::move(recording), lastEvent_);
}

void SuperStepRecorder::beginReplay(const Runtimes& runtimes) {
  waitFor_.clear();
  for (auto* runtime : runtimes) {
    // the replays follow each other on the stream of the recorder; a runtime that has had no work
    // since the latest one has nothing new to wait for
    if (!runtime->defersWaitFor(lastEvent_)) {
      waitFor_.push_back(runtime->eventRecord());
    }
  }
}

void SuperStepRecorder::replay(const Key& key,
                               const std::vector<AbstractTimeCluster*>& clusters,
                               const Runtimes& runtimes,
                               double* deviceTime) {
  launch(graphs_.at(key), clusters, runtimes, deviceTime);
  used_.insert(key);
}

void SuperStepRecorder::forgetUnused() {
  const bool unused = std::any_of(graphs_.begin(), graphs_.end(), [this](const auto& graph) {
    return used_.count(graph.first) == 0;
  });
  if (unused) {
    // all recordings run on the stream of the recorder
    deviceInstance().api().syncStreamWithHost(stream_);
    for (auto graph = graphs_.begin(); graph != graphs_.end();) {
      graph = used_.count(graph->first) == 0 ? graphs_.erase(graph) : std::next(graph);
    }
  }
  used_.clear();
}

void SuperStepRecorder::launch(const device::DeviceGraphHandle& graph,
                               const std::vector<AbstractTimeCluster*>& clusters,
                               const Runtimes& runtimes,
                               double* deviceTime) {
  for (const auto& event : waitFor_) {
    runtime_.eventSync(event);
  }
  waitFor_.clear();
  // the time from where the replay can start
  const bool timed = timer_ != nullptr && deviceTime != nullptr;
  const auto interval = timed ? timer_->begin(stream_) : DeviceTimer::Interval{};
  runtime_.launchGraph(graph);
  if (timed) {
    timer_->end(interval, stream_, deviceTime);
  }

  // everything enqueued from now on comes after the replayed work; the runtimes wait for it only
  // once they get work, and not at all if the next replay comes first
  lastEvent_ = runtime_.eventRecord();
  for (auto* runtime : runtimes) {
    runtime->deferWait(lastEvent_);
  }
  for (auto* cluster : clusters) {
    cluster->publishEvent(lastEvent_);
  }
}

#else

SuperStepRecorder::SuperStepRecorder(bool /*branches*/) {}

void SuperStepRecorder::dispose() {}

bool SuperStepRecorder::available() { return false; }

bool SuperStepRecorder::has(const Key& /*key*/) const { return false; }

ActorEvent SuperStepRecorder::lastEvent() const { return {}; }

void SuperStepRecorder::beginRecording(const std::vector<AbstractTimeCluster*>& /*clusters*/,
                                       const Runtimes& /*runtimes*/) {}

void SuperStepRecorder::endRecording(const Key& /*key*/) {}

void SuperStepRecorder::endRecordingAndLaunch(const std::vector<AbstractTimeCluster*>& /*clusters*/,
                                              const Runtimes& /*runtimes*/,
                                              double* /*deviceTime*/) {}

void SuperStepRecorder::beginReplay(const Runtimes& /*runtimes*/) {}

void SuperStepRecorder::replay(const Key& /*key*/,
                               const std::vector<AbstractTimeCluster*>& /*clusters*/,
                               const Runtimes& /*runtimes*/,
                               double* /*deviceTime*/) {}

void SuperStepRecorder::forgetUnused() {}

#endif

} // namespace seissol::solver

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DeviceLoopTimer.h"

#include "Monitoring/LoopStatistics.h"

#ifdef ACL_DEVICE
#include "Parallel/Runtime/Stream.h"

#include <Device/AbstractAPI.h>
#include <Device/device.h>
#endif

namespace seissol {

namespace {
#ifdef ACL_DEVICE
device::AbstractAPI& api() { return device::DeviceInstance::instance().api(); }
#endif
} // namespace

DeviceLoopTimer::~DeviceLoopTimer() { release(); }

void DeviceLoopTimer::release() {
#ifdef ACL_DEVICE
  collectAll();
  for (const auto& [start, end] : free_) {
    api().destroyEvent(start);
    api().destroyEvent(end);
  }
  free_.clear();
#endif
}

DeviceLoopTimer::Loop DeviceLoopTimer::begin([[maybe_unused]] void* stream) {
  Loop loop;
  if constexpr (!available()) {
    loop.hostBegin = LoopStatistics::begin();
    return loop;
  }
#ifdef ACL_DEVICE
  if (!enabled()) {
    loop.timed = false;
    return loop;
  }
  if (parallel::runtime::recordingOuterGraph()) {
    // the events would become part of the recording, and get recorded anew with each replay
    loop.timed = false;
    return loop;
  }
  if (free_.empty()) {
    loop.events = {api().createEvent(true), api().createEvent(true)};
  } else {
    loop.events = free_.back();
    free_.pop_back();
  }
  api().recordEventOnStream(loop.events.first, stream);
#endif
  return loop;
}

void DeviceLoopTimer::end(const Loop& loop,
                          [[maybe_unused]] void* stream,
                          unsigned region,
                          unsigned numIterations,
                          unsigned subRegion) {
  if (!loop.timed) {
    return;
  }
  if constexpr (!available()) {
    statistics_->end(region, numIterations, subRegion, loop.hostBegin);
    return;
  }
#ifdef ACL_DEVICE
  api().recordEventOnStream(loop.events.second, stream);
  pending_.push_back({loop.events, region, numIterations, subRegion});
#endif
}

void DeviceLoopTimer::collect() {
#ifdef ACL_DEVICE
  if (parallel::runtime::recordingOuterGraph()) {
    // querying an event is not allowed while the streams get captured
    return;
  }
  while (!pending_.empty() && api().isEventCompleted(pending_.front().events.second)) {
    add(pending_.front());
    pending_.pop_front();
  }
#endif
}

void DeviceLoopTimer::collectAll() {
#ifdef ACL_DEVICE
  while (!pending_.empty()) {
    api().syncEventWithHost(pending_.front().events.second);
    add(pending_.front());
    pending_.pop_front();
  }
#endif
}

void DeviceLoopTimer::add([[maybe_unused]] const Pending& loop) {
#ifdef ACL_DEVICE
  const auto time = api().timespanEvents(loop.events.first, loop.events.second);
  statistics_->addTime(loop.region, loop.numIterations, loop.subRegion, time);
  free_.push_back(loop.events);
#endif
}

} // namespace seissol

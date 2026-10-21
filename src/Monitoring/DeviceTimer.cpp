// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DeviceTimer.h"

#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Trace.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

#ifdef ACL_DEVICE
#include "Monitoring/DeviceTimeline.h"
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

DeviceTimer::~DeviceTimer() { release(); }

void DeviceTimer::release() {
#ifdef ACL_DEVICE
  collectAll();
  for (const auto& [start, end] : free_) {
    api().destroyEvent(start);
    api().destroyEvent(end);
  }
  free_.clear();
#endif
}

DeviceTimer::Interval DeviceTimer::begin([[maybe_unused]] void* stream) {
  Interval interval;
  if constexpr (!available()) {
    interval.hostBegin = LoopStatistics::begin();
    return interval;
  }
#ifdef ACL_DEVICE
  if (!enabled()) {
    interval.timed = false;
    return interval;
  }
  if (parallel::runtime::recordingOuterGraph()) {
    // the events would become part of the recording, and get recorded anew with each replay
    interval.timed = false;
    return interval;
  }
  if (free_.empty()) {
    interval.events = {api().createEvent(true), api().createEvent(true)};
  } else {
    interval.events = free_.back();
    free_.pop_back();
  }
  api().recordEventOnStream(interval.events.first, stream);
#endif
  return interval;
}

void DeviceTimer::end(const Interval& interval,
                      void* stream,
                      unsigned region,
                      unsigned numIterations,
                      unsigned subRegion) {
  if (!interval.timed) {
    return;
  }
  if constexpr (!available()) {
    statistics_->end(region, numIterations, subRegion, interval.hostBegin);
    return;
  }
  push(interval, stream, {interval.events, stream, region, numIterations, subRegion, nullptr});
}

void DeviceTimer::end(const Interval& interval, void* stream, double* total) {
  if constexpr (!available()) {
    // the time of the host would be that of the launches
    return;
  }
  if (!interval.timed) {
    return;
  }
  push(interval, stream, {interval.events, stream, 0, 0, 0, total});
}

void DeviceTimer::push([[maybe_unused]] const Interval& interval,
                       [[maybe_unused]] void* stream,
                       [[maybe_unused]] Pending pending) {
#ifdef ACL_DEVICE
  api().recordEventOnStream(interval.events.second, stream);
  pending_.push_back(pending);
#endif
}

void DeviceTimer::collect() {
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

void DeviceTimer::collectAll() {
#ifdef ACL_DEVICE
  while (!pending_.empty()) {
    api().syncEventWithHost(pending_.front().events.second);
    add(pending_.front());
    pending_.pop_front();
  }
#endif
}

void DeviceTimer::add([[maybe_unused]] const Pending& pending) {
#ifdef ACL_DEVICE
  const auto& [start, end] = pending.events;
  const auto time = api().timespanEvents(start, end);
  if (pending.total != nullptr) {
    *pending.total += time;
  } else {
    statistics_->addTime(pending.region, pending.numIterations, pending.subRegion, time);
  }
  if (timeline_ != nullptr) {
    const bool traced = monitoring::Trace::recording();
    timeline_->add(
        kind_, start, end, traced ? traceName(pending) : 0, traced ? traceTrack(pending) : 0);
  }
  free_.push_back(pending.events);
#endif
}

std::uint32_t DeviceTimer::traceName(const Pending& pending) {
  constexpr auto Unnamed = std::numeric_limits<std::uint32_t>::max();
  const bool loop = pending.total == nullptr && statistics_ != nullptr;
  // the regions of the loops, and then the other work
  const std::size_t index = loop ? pending.region : 0;
  auto& names = loop ? traceRegions_ : traceOther_;
  if (names.size() <= index) {
    names.resize(index + 1, Unnamed);
  }
  if (names[index] == Unnamed) {
    names[index] =
        monitoring::Trace::name(loop ? statistics_->regionName(pending.region) : traceWork_);
  }
  return names[index];
}

std::uint32_t DeviceTimer::traceTrack(const Pending& pending) {
  for (const auto& [stream, track] : traceStreams_) {
    if (stream == pending.stream) {
      return track;
    }
  }
  // each stream runs its work in order, so its spans follow each other
  const auto track = monitoring::Trace::track(
      pending.stream, "device: " + (traceTrack_ ? traceTrack_() : std::string("stream")));
  traceStreams_.emplace_back(pending.stream, track);
  return track;
}

} // namespace seissol

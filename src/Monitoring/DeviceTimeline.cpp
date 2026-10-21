// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DeviceTimeline.h"

#include "Monitoring/DeviceTimer.h"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/AbstractAPI.h>
#include <Device/device.h>
#endif

namespace seissol {

namespace {
#ifdef ACL_DEVICE
device::AbstractAPI& api() { return device::DeviceInstance::instance().api(); }
#endif
} // namespace

DeviceTimeline::~DeviceTimeline() { release(); }

bool DeviceTimeline::active() { return DeviceTimer::available() && DeviceTimer::enabled(); }

void DeviceTimeline::beginEpoch() {
  for (auto& intervals : intervals_) {
    intervals.clear();
  }
  if (!active()) {
    return;
  }
#ifdef ACL_DEVICE
  if (stream_ == nullptr) {
    // a stream of its own, idle whenever an event gets recorded on it: the event completes right
    // away
    stream_ = api().createStream();
    reference_ = api().createEvent(true);
    end_ = api().createEvent(true);
  }
  api().recordEventOnStream(reference_, stream_);
  inEpoch_ = true;
#endif
}

double DeviceTimeline::sinceBegin([[maybe_unused]] void* event) const {
#ifdef ACL_DEVICE
  return api().timespanEvents(reference_, event);
#else
  return 0;
#endif
}

void DeviceTimeline::add(Kind kind, void* start, void* end) {
  if (!inEpoch_) {
    return;
  }
  intervals_[static_cast<std::size_t>(kind)].emplace_back(sinceBegin(start), sinceBegin(end));
}

DeviceTimeline::Activity DeviceTimeline::endEpoch() {
  Activity activity;
  if (!inEpoch_) {
    return activity;
  }
  inEpoch_ = false;
#ifdef ACL_DEVICE
  api().recordEventOnStream(end_, stream_);
  api().syncEventWithHost(end_);
  activity.span = api().timespanEvents(reference_, end_);
#endif
  std::vector<std::pair<double, double>> all;
  for (std::size_t kind = 0; kind < KindCount; ++kind) {
    activity.busyBy[kind] = unionLength(intervals_[kind]);
    all.insert(all.end(), intervals_[kind].begin(), intervals_[kind].end());
  }
  activity.busy = unionLength(all);
  // an exchange may wait for another process on the device (e.g. NCCL), which keeps the device
  // busy without computing
  const auto& exchanges = intervals_[static_cast<std::size_t>(Kind::Exchanges)];
  std::vector<std::pair<double, double>> computing;
  computing.reserve(all.size() - exchanges.size());
  for (const auto kind : {Kind::Loops, Kind::SuperSteps}) {
    const auto& intervals = intervals_[static_cast<std::size_t>(kind)];
    computing.insert(computing.end(), intervals.begin(), intervals.end());
  }
  activity.computing = unionLength(std::move(computing));
  return activity;
}

void DeviceTimeline::release() {
#ifdef ACL_DEVICE
  if (stream_ != nullptr) {
    api().destroyEvent(reference_);
    api().destroyEvent(end_);
    api().destroyGenericStream(stream_);
    stream_ = nullptr;
    reference_ = nullptr;
    end_ = nullptr;
  }
#endif
  inEpoch_ = false;
}

double DeviceTimeline::unionLength(std::vector<std::pair<double, double>> intervals) {
  std::sort(intervals.begin(), intervals.end());
  double length = 0;
  double coveredUntil = 0;
  bool first = true;
  for (const auto& [from, to] : intervals) {
    if (first || from > coveredUntil) {
      length += to - from;
      coveredUntil = to;
      first = false;
    } else if (to > coveredUntil) {
      length += to - coveredUntil;
      coveredUntil = to;
    }
  }
  return length;
}

} // namespace seissol

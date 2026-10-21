// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_DEVICETIMER_H_
#define SEISSOL_SRC_MONITORING_DEVICETIMER_H_

#include "Monitoring/DeviceTimeline.h"
#include "Monitoring/LoopStatistics.h"

#include <deque>
#include <time.h>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include "Alignment.h"
#include "Common/Typedefs.h"
#endif

namespace seissol {

/**
 * Times the work enqueued on a device stream: the loops of a cluster for the loop statistics, or
 * e.g. the replays of super-timesteps for a total. A pair of timing events goes around each piece
 * of work on the stream, and its time counts once the device is past it, also in the timeline of
 * the epoch (see DeviceTimeline). Taking the time never waits for the device; collect() adds the
 * work that is done.
 *
 * The time of the host, which only enqueues the work, would be the time of the launches. The SYCL
 * backend times its events only when the device library is built for profiling, so the loops keep
 * the time of the host there. Nothing gets timed within a recording of super-timesteps; the work
 * replayed from it counts as the replay.
 */
class DeviceTimer {
  public:
  DeviceTimer(LoopStatistics* statistics, DeviceTimeline* timeline, DeviceTimeline::Kind kind)
      : statistics_(statistics), timeline_(timeline), kind_(kind) {}
  ~DeviceTimer();

  DeviceTimer(const DeviceTimer&) = delete;
  DeviceTimer(DeviceTimer&&) = delete;
  DeviceTimer& operator=(const DeviceTimer&) = delete;
  DeviceTimer& operator=(DeviceTimer&&) = delete;

  /// whether this build times the work on the device
  static constexpr bool available();

  /**
   * Whether the work on the device gets timed (SEISSOL_DEVICE_LOOP_STATISTICS); otherwise, it does
   * not count at all. Each timed piece of work costs two records of events.
   */
  static bool& enabled() {
    static bool enabled = true;
    return enabled;
  }

  /**
   * A piece of work being timed: by a pair of events, by the time of the host, or not at all
   * (within a recording of super-timesteps).
   */
  struct Interval {
    std::pair<void*, void*> events{nullptr, nullptr};
    timespec hostBegin{};
    bool timed{true};
  };

  /**
   * Starts timing work that gets enqueued on `stream` from now on.
   */
  [[nodiscard]] Interval begin(void* stream);

  /**
   * Ends timing a loop, enqueued on `stream` up to now; its time counts for the loop statistics
   * once the device is past it.
   */
  void end(const Interval& interval,
           void* stream,
           unsigned region,
           unsigned numIterations,
           unsigned subRegion);

  /**
   * Ends timing work, enqueued on `stream` up to now; its time gets added to `total` once the
   * device is past it.
   */
  void end(const Interval& interval, void* stream, double* total);

  /**
   * Adds the times of the work the device is past, in the order it was enqueued; without waiting,
   * and not while a recording of super-timesteps captures the streams.
   */
  void collect();

  /**
   * Waits until the device is past all work, and adds its times; e.g. at the end of an epoch.
   */
  void collectAll();

  /**
   * Adds the times of all work, and gives the events back; before the device goes away.
   */
  void release();

  private:
  struct Pending {
    std::pair<void*, void*> events;
    unsigned region{0};
    unsigned numIterations{0};
    unsigned subRegion{0};
    // where the time goes instead of the loop statistics
    double* total{nullptr};
  };

  /// records the end of the work, which counts as `pending` then
  void push(const Interval& interval, void* stream, Pending pending);

  /// adds the time of the oldest work
  void add(const Pending& pending);

  LoopStatistics* statistics_;
  // (only with a device)
  [[maybe_unused]] DeviceTimeline* timeline_;
  [[maybe_unused]] DeviceTimeline::Kind kind_;
  std::deque<Pending> pending_;
  std::vector<std::pair<void*, void*>> free_;
};

constexpr bool DeviceTimer::available() {
#ifdef ACL_DEVICE
  return Backend == DeviceBackend::Cuda || Backend == DeviceBackend::Hip;
#else
  return false;
#endif
}

} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_DEVICETIMER_H_

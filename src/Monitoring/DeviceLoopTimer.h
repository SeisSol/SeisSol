// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_DEVICELOOPTIMER_H_
#define SEISSOL_SRC_MONITORING_DEVICELOOPTIMER_H_

#include "Alignment.h"
#include "Common/Typedefs.h"
#include "Monitoring/LoopStatistics.h"

#include <deque>
#include <time.h>
#include <utility>
#include <vector>

namespace seissol {

/**
 * Times the loops that a cluster enqueues on its device stream, for the loop statistics: a pair of
 * timing events goes around each loop on the stream, and its time counts once the device is past
 * it. Taking the time never waits for the device; collect() adds the loops that are done.
 *
 * The time of the host, which only enqueues the loops, would be the time of the launches. The SYCL
 * backend times its events only when the device library is built for profiling, so the loops keep
 * the time of the host there. The steps replayed from a recording of super-timesteps do not go
 * through the loops; their time does not count.
 */
class DeviceLoopTimer {
  public:
  explicit DeviceLoopTimer(LoopStatistics* statistics) : statistics_(statistics) {}
  ~DeviceLoopTimer();

  DeviceLoopTimer(const DeviceLoopTimer&) = delete;
  DeviceLoopTimer(DeviceLoopTimer&&) = delete;
  DeviceLoopTimer& operator=(const DeviceLoopTimer&) = delete;
  DeviceLoopTimer& operator=(DeviceLoopTimer&&) = delete;

  /// whether this build times the loops on the device
  static constexpr bool available();

  /**
   * Whether the loops on the device get timed (SEISSOL_DEVICE_LOOP_STATISTICS); otherwise, they do
   * not count at all. Each timed loop costs two records of events.
   */
  static bool& enabled() {
    static bool enabled = true;
    return enabled;
  }

  /**
   * A loop being timed: by a pair of events, by the time of the host, or not at all (within a
   * recording of super-timesteps).
   */
  struct Loop {
    std::pair<void*, void*> events{nullptr, nullptr};
    timespec hostBegin{};
    bool timed{true};
  };

  /**
   * Starts timing a loop that gets enqueued on `stream` from now on.
   */
  [[nodiscard]] Loop begin(void* stream);

  /**
   * Ends timing the loop, enqueued on `stream` up to now. Its time counts once the device is past
   * it.
   */
  void end(
      const Loop& loop, void* stream, unsigned region, unsigned numIterations, unsigned subRegion);

  /**
   * Adds the times of the loops the device is past, in the order they were enqueued; without
   * waiting, and not while a recording of super-timesteps captures the streams.
   */
  void collect();

  /**
   * Waits until the device is past all loops, and adds their times; e.g. at the end of an epoch.
   */
  void collectAll();

  /**
   * Adds the times of all loops, and gives the events back; before the device goes away.
   */
  void release();

  private:
  struct Pending {
    std::pair<void*, void*> events;
    unsigned region;
    unsigned numIterations;
    unsigned subRegion;
  };

  /// adds the time of the oldest loop
  void add(const Pending& loop);

  LoopStatistics* statistics_;
  std::deque<Pending> pending_;
  std::vector<std::pair<void*, void*>> free_;
};

constexpr bool DeviceLoopTimer::available() {
#ifdef ACL_DEVICE
  return Backend == DeviceBackend::Cuda || Backend == DeviceBackend::Hip;
#else
  return false;
#endif
}

} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_DEVICELOOPTIMER_H_

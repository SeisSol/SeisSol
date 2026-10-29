// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_DEVICETIMELINE_H_
#define SEISSOL_SRC_MONITORING_DEVICETIMELINE_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace seissol {

/**
 * What the device did in an epoch, from the timing events around the work enqueued on it (see
 * DeviceTimer): the intervals in which it worked, by kind, as times relative to an event at the
 * start of the epoch. Their union tells how long the device was busy in the epoch, and with that,
 * how long it waited: e.g. for the host to enqueue the work, for other processes, or at the
 * synchronization point.
 *
 * The intervals of a kind may overlap (e.g. the loops of clusters that run concurrently); the
 * union counts each moment once.
 */
class DeviceTimeline {
  public:
  enum class Kind : std::uint8_t {
    /// the loops of the clusters
    Loops,
    /// the replays of recorded super-timesteps
    SuperSteps,
    /// the exchanges of the halo on device streams
    Exchanges
  };
  static constexpr std::size_t KindCount = 3;

  /**
   * The device time of an epoch.
   */
  struct Activity {
    /// from the start of the epoch until the device was done with all its work
    double span{0};
    /// the time in which the device worked on anything
    double busy{0};
    /// the time in which it computed: ran loops or replays of super-timesteps
    double computing{0};
    /// the time in which the device worked on each kind
    std::array<double, KindCount> busyBy{};

    [[nodiscard]] double busyWith(Kind kind) const {
      return busyBy[static_cast<std::size_t>(kind)];
    }

    /// adds the times of another epoch
    Activity& operator+=(const Activity& other) {
      span += other.span;
      busy += other.busy;
      computing += other.computing;
      for (std::size_t kind = 0; kind < KindCount; ++kind) {
        busyBy[kind] += other.busyBy[kind];
      }
      return *this;
    }
  };

  DeviceTimeline() = default;
  ~DeviceTimeline();

  DeviceTimeline(const DeviceTimeline&) = delete;
  DeviceTimeline(DeviceTimeline&&) = delete;
  DeviceTimeline& operator=(const DeviceTimeline&) = delete;
  DeviceTimeline& operator=(DeviceTimeline&&) = delete;

  /// whether the work on the device gets timed (see DeviceTimer::available() and enabled())
  [[nodiscard]] static bool active();

  /**
   * With a trace (see monitoring::Trace) and SEISSOL_TRACE_KERNELS=1: lets the device put each of
   * its kernels and waits into it, on a track for each of its queues, where it can (see
   * device::AbstractAPI::traceDevice(); HIP with own queues). stopKernelTrace() waits until all of
   * them are in; call it before the trace gets written.
   */
  static void startKernelTrace(bool requested);
  static void stopKernelTrace();

  /**
   * Starts an epoch; once the device is done with all work before it.
   */
  void beginEpoch();

  /**
   * Adds the interval between two events of the current epoch, which have completed; while a trace
   * records (see monitoring::Trace), also as a span with the given name on the given track.
   */
  void add(Kind kind, void* start, void* end, std::uint32_t traceName, std::uint32_t traceTrack);

  /**
   * Ends the epoch, once the device is done with all of its work, and all intervals have been
   * added.
   */
  [[nodiscard]] Activity endEpoch();

  /// gives back the events and the stream; before the device goes away
  void release();

  /// the length of the union of the intervals
  [[nodiscard]] static double unionLength(std::vector<std::pair<double, double>> intervals);

  private:
  /// the time of a completed event of the epoch, relative to its start
  [[nodiscard]] double sinceBegin(void* event) const;

  // (only with a device)
  [[maybe_unused]] void* stream_{nullptr};
  [[maybe_unused]] void* reference_{nullptr};
  [[maybe_unused]] void* end_{nullptr};
  bool inEpoch_{false};
  // in a trace: the time of the start of the epoch on the host, in microseconds since its origin
  [[maybe_unused]] double traceBegin_{0};
  std::array<std::vector<std::pair<double, double>>, KindCount> intervals_;
};

} // namespace seissol

#endif // SEISSOL_SRC_MONITORING_DEVICETIMELINE_H_

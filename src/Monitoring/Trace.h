// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_TRACE_H_
#define SEISSOL_SRC_MONITORING_TRACE_H_

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mpi.h>
#include <string>
#include <string_view>
#include <utils/env.h>

namespace seissol::monitoring {

/**
 * A timeline of what the processes do, for viewers such as Perfetto (https://ui.perfetto.dev) or
 * chrome://tracing: spans of time on tracks, written in the trace event format of Chrome, as one
 * JSON file for all processes. Each process has its tracks: the threads of the host, the streams
 * of the device, and others, e.g. for the halo exchanges.
 *
 * Recording costs little, and nothing but a check while the trace is off: the spans go into
 * buffers of the threads that record them, and into the file at the end. The trace records only
 * within its window of time, and at most a given number of spans per process.
 */
class Trace {
  public:
  struct Settings {
    bool enabled{false};
    /// the window, in seconds since the start of the time stepping
    double start{0};
    double duration{std::numeric_limits<double>::infinity()};
    /// the most spans a process records
    std::size_t maxSpans{1000000};
  };

  /// SEISSOL_TRACE, SEISSOL_TRACE_START and SEISSOL_TRACE_DURATION (in seconds),
  /// SEISSOL_TRACE_SPANS
  static Settings settingsFromEnv(utils::Env& env);

  /**
   * Starts recording, with the same origin of time on all processes (up to the skew of a
   * barrier); nothing to do if the trace is off. Collective.
   */
  static void start(const Settings& settings, MPI_Comm comm);

  /// whether the trace records
  [[nodiscard]] static bool recording() { return recordingFlag().load(std::memory_order_acquire); }

  /// a point in time, in microseconds since the origin
  [[nodiscard]] static double time(std::chrono::steady_clock::time_point point);
  [[nodiscard]] static double now() { return time(std::chrono::steady_clock::now()); }

  /// the id of a name of spans; the same for the same name
  [[nodiscard]] static std::uint32_t name(std::string_view name);

  /**
   * A track of its own, e.g. for a stream of the device; the same for the same key, and named by
   * the first name given with it.
   */
  [[nodiscard]] static std::uint32_t track(const void* key, std::string_view name);

  /// the track of the calling thread
  [[nodiscard]] static std::uint32_t threadTrack();

  /**
   * Records a span on a track, with its begin and end in microseconds since the origin; from any
   * thread. A track shows its spans nested: two spans of a track must not overlap partially.
   */
  static void span(std::uint32_t name, std::uint32_t track, double begin, double end);

  /**
   * Stops recording, and writes the trace of all processes into `path`; nothing to do if the
   * trace is off. Collective.
   */
  static void write(const std::string& path, MPI_Comm comm);

  private:
  static std::atomic<bool>& recordingFlag();
};

/**
 * Records a span on the track of the calling thread, from its construction to its destruction.
 */
class TraceScope {
  public:
  explicit TraceScope(std::uint32_t name)
      : name_(name), begin_(Trace::recording() ? Trace::now() : -1) {}
  ~TraceScope() {
    if (begin_ >= 0 && Trace::recording()) {
      Trace::span(name_, Trace::threadTrack(), begin_, Trace::now());
    }
  }

  TraceScope(const TraceScope&) = delete;
  TraceScope(TraceScope&&) = delete;
  TraceScope& operator=(const TraceScope&) = delete;
  TraceScope& operator=(TraceScope&&) = delete;

  private:
  std::uint32_t name_;
  double begin_;
};

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_TRACE_H_

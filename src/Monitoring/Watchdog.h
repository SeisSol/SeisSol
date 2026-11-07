// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_WATCHDOG_H_
#define SEISSOL_SRC_MONITORING_WATCHDOG_H_

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <utils/env.h>
#include <vector>

namespace seissol::monitoring {

/**
 * Watches the progress of the simulation from a thread of its own, also where the main thread
 * cannot report it any more, e.g. when it hangs in a call to the device or to a communication
 * library.
 *
 * Probes are values that change with the progress (e.g. the step counts of the clusters, or the
 * times the device writes back for them), read without synchronization with whoever changes
 * them. While the simulation takes its steps, none of them changing for the timeout is a stall; in
 * the other phases (e.g. the outputs), the phase lasting longer than the timeout is one. Checks run
 * at every period, e.g. for the asynchronous errors of a library.
 *
 * At a stall, the watchdog reports the phase and the probes, and asks for the reports of those who
 * know more (see addReport()). A period later, it aborts a stall in the time stepping, if
 * configured to; a phase elsewhere that lasts long (e.g. a large checkpoint, or the last writes of
 * the outputs) it only reports, as writing much data may well take that long.
 */
class Watchdog {
  public:
  struct Settings {
    bool enabled{true};
    std::chrono::seconds timeout{900};
    std::chrono::seconds period{10};
    bool abort{true};
  };

  /// SEISSOL_WATCHDOG, SEISSOL_WATCHDOG_TIMEOUT (in seconds), SEISSOL_WATCHDOG_ABORT
  static Settings settingsFromEnv(utils::Env& env);

  enum class Phase : std::uint8_t {
    //! not watched (e.g. before the simulation)
    Idle,
    //! the time stepping: stalled if no probe changes for the timeout
    Stepping,
    //! anything else: stalled if it lasts longer than the timeout; reported, but never aborted
    Other
  };

  /// What a look at the progress has found.
  enum class Finding : std::uint8_t {
    //! progress, or nothing watched
    Progress,
    //! a stall, reported once per phase
    Stall,
    //! a stall in the time stepping, a period after its report: the simulation is to be aborted
    Abort
  };

  explicit Watchdog(const Settings& settings);
  ~Watchdog();

  Watchdog(const Watchdog&) = delete;
  Watchdog& operator=(const Watchdog&) = delete;
  Watchdog(Watchdog&&) = delete;
  Watchdog& operator=(Watchdog&&) = delete;

  /// Adds a probe; before start().
  void addProbe(const std::string& name, std::function<double()> value);

  /// Adds a check that runs at every period; before start().
  void addCheck(std::function<void()> check);

  /// Adds what asks for a report at a stall (e.g. of a thread that may still run); before start().
  void addReport(std::function<void()> report);

  void start();
  void stop();

  /// What the main thread does now; `label` names it in the reports, and has to stay valid.
  void setPhase(Phase phase, const char* label);

  [[nodiscard]] const Settings& settings() const { return settings_; }

  using Clock = std::chrono::steady_clock;

  /**
   * One look at the progress at the time `now`, as the thread of the watchdog takes one every
   * period: runs the checks, reads the probes, and reports a stall. Whether to abort, the thread
   * decides by what it finds and the settings.
   */
  Finding watch(Clock::time_point now);

  private:
  struct Probe {
    std::string name;
    std::function<double()> value;
    double last{0};
    Clock::time_point lastChange;
  };

  void run();

  Settings settings_;
  std::vector<Probe> probes_;
  std::vector<std::function<void()>> checks_;
  std::vector<std::function<void()>> reports_;

  std::atomic<Phase> phase_{Phase::Idle};
  std::atomic<const char*> label_{"idle"};
  std::atomic<std::int64_t> phaseStart_{0};
  std::atomic<std::uint64_t> phaseCount_{0};

  // the phase a stall was reported for, so that it gets reported once
  std::uint64_t reportedPhase_{0};
  bool reported_{false};

  std::mutex mutex_;
  std::condition_variable wakeup_;
  bool stop_{false};
  std::thread thread_;
};

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_WATCHDOG_H_

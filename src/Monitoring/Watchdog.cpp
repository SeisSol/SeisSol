// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Watchdog.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <utils/env.h>
#include <utils/logger.h>

namespace seissol::monitoring {

Watchdog::Settings Watchdog::settingsFromEnv(utils::Env& env) {
  Settings settings;
  settings.enabled = env.get<bool>("WATCHDOG", true);
  const auto timeout = env.get<double>("WATCHDOG_TIMEOUT", 900.0);
  settings.timeout =
      std::chrono::seconds(std::max<std::int64_t>(1, static_cast<std::int64_t>(timeout)));
  // at least ten looks per timeout, and one every ten seconds at most
  settings.period =
      std::clamp(settings.timeout / 10, std::chrono::seconds(1), std::chrono::seconds(10));
  settings.abort = env.get<bool>("WATCHDOG_ABORT", true);
  return settings;
}

Watchdog::Watchdog(const Settings& settings) : settings_(settings) {}

Watchdog::~Watchdog() { stop(); }

void Watchdog::addProbe(const std::string& name, std::function<double()> value) {
  assert(!thread_.joinable());
  probes_.push_back(Probe{name, std::move(value), 0, Clock::now()});
}

void Watchdog::addCheck(std::function<void()> check) {
  assert(!thread_.joinable());
  checks_.push_back(std::move(check));
}

void Watchdog::addReport(std::function<void()> report) {
  assert(!thread_.joinable());
  reports_.push_back(std::move(report));
}

void Watchdog::start() {
  if (!settings_.enabled || thread_.joinable()) {
    return;
  }
  const auto now = Clock::now();
  for (auto& probe : probes_) {
    probe.last = probe.value();
    probe.lastChange = now;
  }
  stop_ = false;
  thread_ = std::thread([this] { run(); });
}

void Watchdog::stop() {
  if (!thread_.joinable()) {
    return;
  }
  {
    const std::scoped_lock lock(mutex_);
    stop_ = true;
  }
  wakeup_.notify_all();
  thread_.join();
}

void Watchdog::setPhase(Phase phase, const char* label) {
  label_ = label;
  phaseStart_ = Clock::now().time_since_epoch().count();
  phase_ = phase;
  ++phaseCount_;
}

void Watchdog::run() {
  std::unique_lock lock(mutex_);
  while (!wakeup_.wait_for(lock, settings_.period, [this] { return stop_; })) {
    lock.unlock();
    if (watch(Clock::now()) == Finding::Abort && settings_.abort) {
      logError() << "Watchdog: aborting the stalled simulation (SEISSOL_WATCHDOG_ABORT=0 keeps it"
                 << "waiting).";
    }
    lock.lock();
  }
}

Watchdog::Finding Watchdog::watch(Clock::time_point now) {
  for (const auto& check : checks_) {
    check();
  }

  const auto phase = phase_.load();
  const auto phaseCount = phaseCount_.load();
  const auto phaseStart = Clock::time_point(Clock::duration(phaseStart_.load()));

  // in the time stepping, the last progress of any probe counts; elsewhere, the phase itself
  auto since = phaseStart;
  for (auto& probe : probes_) {
    const auto value = probe.value();
    if (value != probe.last) {
      probe.last = value;
      probe.lastChange = now;
    }
    if (phase == Phase::Stepping) {
      since = std::max(since, probe.lastChange);
    }
  }

  const bool stalled = phase != Phase::Idle && now - since > settings_.timeout;
  if (!stalled) {
    reported_ = false;
    return Finding::Progress;
  }

  const auto seconds = [](auto duration) {
    return std::chrono::duration_cast<std::chrono::seconds>(duration).count();
  };
  if (!reported_ || reportedPhase_ != phaseCount) {
    reported_ = true;
    reportedPhase_ = phaseCount;
    if (phase == Phase::Stepping) {
      logWarning(true) << "Watchdog: no progress for" << seconds(now - since) << "s in"
                       << label_.load() << "(timeout" << seconds(settings_.timeout) << "s).";
    } else {
      logWarning(true) << "Watchdog:" << label_.load() << "has taken" << seconds(now - since)
                       << "s so far, longer than the timeout of" << seconds(settings_.timeout)
                       << "s; outside of the time stepping, the watchdog only reports that.";
    }
    for (const auto& probe : probes_) {
      logWarning(true) << "Watchdog:" << probe.name << "at" << probe.last << ", last changed"
                       << seconds(now - probe.lastChange) << "s ago";
    }
    // the threads that still run add what only they can read
    for (const auto& report : reports_) {
      report();
    }
    return Finding::Stall;
  }

  // a period after the report, for the main thread to add to it; outside of the time stepping, a
  // phase that takes long (e.g. a large checkpoint) is no reason to stop the run
  return phase == Phase::Stepping ? Finding::Abort : Finding::Stall;
}

} // namespace seissol::monitoring

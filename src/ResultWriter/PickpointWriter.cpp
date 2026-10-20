// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PickpointWriter.h"

#include "Modules/Modules.h"

#include <functional>
#include <limits>
#include <optional>
#include <utils/logger.h>

namespace seissol::writer {

void PickpointWriter::captureSyncPoint(double /*currentTime*/) { capture_(); }

void PickpointWriter::processSyncPoint(double /*currentTime*/) {
  logInfo() << "Writing pickpoint data.";

  // currently, we just wrap the DR output function.
  write_();

  logInfo() << "Writing pickpoint data complete.";
}

void PickpointWriter::simulationStart(std::optional<double> /*checkpointTime*/) {}

void PickpointWriter::simulationEnd() { syncPoint(std::numeric_limits<double>::infinity()); }

void PickpointWriter::enable(double interval) {
  if (!enabled_) {
    Modules::registerHook(*this, ModuleHook::SimulationStart);
    Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
    Modules::registerHook(*this, ModuleHook::SimulationEnd);

    setSyncInterval(interval);
  }

  enabled_ = true;
}

void PickpointWriter::setupWriter(const std::function<void()>& capture,
                                  const std::function<void()>& write) {
  capture_ = capture;
  write_ = write;
}

} // namespace seissol::writer

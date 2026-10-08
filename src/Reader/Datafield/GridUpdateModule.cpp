// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Reader/Datafield/GridUpdateModule.h"

#include "Modules/Modules.h"
#include "Reader/Datafield/Grid.h"
#include "utils/logger.h"

#include <optional>
#include <stdexcept>

namespace seissol::reader::datafield {

void GridUpdateModule::simulationStart(std::optional<double> checkpointTime) {
  // On a restart the first query is at the checkpoint time, not at zero, and a
  // window sized for t = 0 would not cover it. Reloading here rather than on
  // the first sample keeps the I/O inside setup, where it is measured, instead
  // of inside the first timestep, where it looks like a solver stall.
  store_->update(checkpointTime.value_or(0.0));
}

void GridUpdateModule::syncPoint(double currentTime) { store_->update(currentTime); }

void GridUpdateModule::enable(double interval) {
  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  // An upper bound on a safe interval, not a request: Modules takes the minimum
  // across every registered module, so a shorter output cadence simply wins and
  // the windows are refreshed more often than strictly needed.
  setSyncInterval(interval);
}

void registerGridUpdateModule() {
  static GridUpdateModule module(sharedGridStore());
  static bool registered = false;
  if (registered) {
    return;
  }

  std::optional<double> interval;
  try {
    interval = sharedGridStore().suggestedSyncInterval();
  } catch (const std::invalid_argument& error) {
    // A window narrower than the interpolation stencil: no interval is safe, so
    // there is nothing to fall back to. The fix is the memory budget or the
    // scheme, both of which the message has to name for the user to act.
    logError() << "A time-dependent grid cannot be advanced safely:" << error.what()
               << "-- raise the window memory budget or use a narrower interpolation in time.";
    return;
  }

  if (!interval.has_value()) {
    logInfo() << "No time-dependent data grids; not registering the grid update module.";
    return;
  }

  module.enable(*interval);
  registered = true;

  logInfo() << "Registered the grid update module with a synchronisation interval of" << *interval
            << "seconds.";
}

} // namespace seissol::reader::datafield

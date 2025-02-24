// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Module.h"

#include "utils/logger.h"
#include <cassert>
#include <cmath>
#include <limits>

namespace seissol {
Module::Module() : lastSyncPoint(-std::numeric_limits<double>::infinity()) {}

Module::~Module() = default;

double Module::potentialSyncPoint(double currentTime, double timeTolerance, bool forceSyncPoint) {
  if (std::abs(currentTime - lastSyncPoint) < timeTolerance) {
    logInfo() << "Ignoring duplicate synchronization point at time" << currentTime
              << "; the last sync point was at " << lastSyncPoint;
  } else if (forceSyncPoint || std::abs(currentTime - nextSyncPoint) < timeTolerance) {
    syncPoint(currentTime);
    lastSyncPoint = currentTime;
    nextSyncPoint += isyncInterval;
  }

  return nextSyncPoint;
}

void Module::setSimulationStartTime(double time) {
  assert(isyncInterval > 0);
  lastSyncPoint = time;

  // take the next expected sync point TODO: forward tolerance
  // (calculated from time point 0)
  nextSyncPoint = 0;
  while (nextSyncPoint - time < 1e-6) {
    nextSyncPoint += isyncInterval;
  }
}

/**
 * Set the synchronization interval for this module
 *
 * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
 */
void Module::setSyncInterval(double interval) {
  if (isyncInterval != 0) {
    logError() << "Synchronization interval is already set";
  }
  isyncInterval = interval;
}
} // namespace seissol

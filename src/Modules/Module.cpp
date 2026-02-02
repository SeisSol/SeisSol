// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Module.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <utils/logger.h>

namespace seissol {
Module::Module() : lastSyncPoint_(-std::numeric_limits<double>::infinity()) {}

Module::~Module() = default;

double Module::potentialSyncPoint(double currentTime, double timeTolerance, bool forceSyncPoint) {
  if (std::abs(currentTime - lastSyncPoint_) < timeTolerance) {
    logDebug() << "Ignoring duplicate synchronization point at time" << currentTime
               << "; the last sync point was at " << lastSyncPoint_;
  } else if (forceSyncPoint || std::abs(currentTime - nextSyncPoint_) < timeTolerance) {
    syncPoint(currentTime);
    lastSyncPoint_ = currentTime;
    nextSyncPoint_ += isyncInterval_;
  }

  return nextSyncPoint_;
}

void Module::setSimulationStartTime(double time) {
  lastSyncPoint_ = time;

  // take the next expected sync point TODO: forward tolerance
  // (calculated from time point 0)
  nextSyncPoint_ = 0;
  while (nextSyncPoint_ - time < 1e-6) {
    const auto currNextSyncPoint = nextSyncPoint_;
    nextSyncPoint_ += isyncInterval_;
    if (currNextSyncPoint == nextSyncPoint_) {
      // no time advancement (i.e. 0 or too small)
      // just jump to the init time
      nextSyncPoint_ = time;
      break;
    }
  }
}

/**
 * Set the synchronization interval for this module
 *
 * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
 */
void Module::setSyncInterval(double interval) {
  if (isyncInterval_ != 0) {
    logError() << "Synchronization interval is already set";
  }
  isyncInterval_ = interval;
}
} // namespace seissol

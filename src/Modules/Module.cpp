// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Module.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <set>
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
    nextSyncPoint_ = nextCandidate(currentTime);
  }

  return nextSyncPoint_;
}

double Module::nextCandidate(double time) const {
  double nextTime = time + syncInterval_;
  if (syncInterval_ > 0) {
    // adjust to a (strict) multiple of the sync interval; just to be sure
    nextTime = std::floor((time + syncInterval_) / syncInterval_) * syncInterval_;
  }

  for (const auto& extra : extraPoints_) {
    if (extra > time) {
      nextTime = std::min(extra, nextTime);
    }
  }

  return nextTime;
}

void Module::setSimulationStartTime(double time) {
  lastSyncPoint_ = time;

  nextSyncPoint_ = nextCandidate(time);
}

/**
 * Set the synchronization interval for this module
 *
 * This is only required for modules that register for {@link SYNCHRONIZATION_POINT}.
 */
void Module::setSyncInterval(double interval) {
  if (syncInterval_ != 0) {
    logError() << "Synchronization interval is already set";
  }
  syncInterval_ = interval;
}

void Module::addExtraSyncPoints(const std::vector<double>& points) {
  std::set<double> allpoints(extraPoints_.begin(), extraPoints_.end());
  allpoints.insert(points.begin(), points.end());
  extraPoints_ = std::vector<double>(allpoints.begin(), allpoints.end());
}

} // namespace seissol

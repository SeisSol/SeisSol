// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "FrictionSolver.h"

#include "Config.h"
#include "Features/DynamicRupture/Misc.h"

#include <cstddef>
#include <utils/logger.h>
#include <vector>

namespace seissol::dr::friction_law {

template <typename Cfg>
FrictionSolver::FrictionTime FrictionSolver::computeDeltaT(const std::vector<double>& timePoints) {
  std::vector<double> deltaT(misc::TimeSteps<Cfg>);

  if (timePoints.size() != deltaT.size()) {
    logError() << "Internal time point count mismatch. Given vs. expected:" << timePoints.size()
               << deltaT.size();
  }

  deltaT[0] = timePoints[0]; // - 0
  for (std::size_t timeIndex = 1; timeIndex < misc::TimeSteps<Cfg>; ++timeIndex) {
    deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
  }

  // add the segment [lastPoint, timestep] to the last point
  deltaT.back() += timePoints[0];

  return {deltaT};
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template FrictionSolver::FrictionTime FrictionSolver::computeDeltaT<Cfg>(                        \
      const std::vector<double>& timePoints);                                                      \
  template class FrictionSolverImpl<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law

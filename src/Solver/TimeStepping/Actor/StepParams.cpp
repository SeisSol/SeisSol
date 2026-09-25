// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "StepParams.h"

#include "Solver/TimeStepping/Actor/ActorState.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace seissol::time_stepping {

StepContext computeStepContext(const ClusterTimes& times,
                               const std::vector<NeighborCluster>& neighbors) {
  StepContext context{};
  context.largestTimeStepSize = times.maxTimeStepSize;
  for (const auto& neighbor : neighbors) {
    context.largestTimeStepSize =
        std::max(context.largestTimeStepSize, neighbor.ct.maxTimeStepSize);
    if (neighbor.ct.timeStepRate > times.timeStepRate) {
      // neighboring clusters differ by at most one cluster index, so all larger neighbors share
      // one update rate
      assert(context.largerNeighborRate == 0 ||
             context.largerNeighborRate == neighbor.ct.timeStepRate);
      context.largerNeighborRate = neighbor.ct.timeStepRate;
    }
  }
  return context;
}

StepParams
    computeStepParams(const ClusterTimes& times, double syncTime, const StepContext& context) {
  StepParams params{};
  params.time = times.correctionTime;
  params.timeStepSize = times.timeStepSize(syncTime);
  params.neighborTimeStepSize = context.largestTimeStepSize;

  if (context.largerNeighborRate > 0) {
    assert(times.timeStepRate > 0);
    assert(context.largerNeighborRate % times.timeStepRate == 0);

    // both clusters restart their step counts at each synchronization point, so the position
    // inside the neighbor's step follows from the own step count alone
    const auto substep =
        (times.stepsSinceLastSync % context.largerNeighborRate) / times.timeStepRate;
    params.subTimeStart = static_cast<double>(substep) * times.maxTimeStepSize;
    params.resetBuffers = substep == 0;
  }

  return params;
}

} // namespace seissol::time_stepping

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_TIMESTEPPINGPLAN_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_TIMESTEPPINGPLAN_H_

#include "Solver/TimeStepping/Actor/ActorState.h"

#include <cstddef>
#include <vector>

namespace seissol::solver {

/**
 * What the plan needs to know about a cluster.
 */
struct PlannedCluster {
  long timeStepRate{1};
  long stepsUntilSync{0};
  DataReadiness dataReadiness{DataReadiness::AfterPrediction};
  ActorPriority priority{ActorPriority::Low};
};

/**
 * One action of one cluster in the plan.
 */
struct PlannedAction {
  std::size_t cluster;
  ActorAction action;
  long step;
};

/**
 * Fixes in advance the order in which the clusters of a process take their steps between two
 * synchronization points.
 *
 * Every action is placed at a point in logical time, counted in steps of the smallest time cluster
 * since the synchronization point: the correction of a cell cluster at the end of its step, every
 * other action at the start of its step. At the same point, the corrections come first, since they
 * complete the steps before it; then the predictions, which start the steps from there on; then the
 * face work, which needs the predictions on both sides. Everything a cluster waits for thereby
 * comes earlier in the plan. Within these rules, clusters of high priority go first.
 *
 * The plan leaves out the ghost clusters: their progress depends on other processes.
 */
std::vector<PlannedAction> planTimeSteps(const std::vector<PlannedCluster>& clusters);

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_PLAN_TIMESTEPPINGPLAN_H_

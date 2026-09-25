// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_STEPPARAMS_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_STEPPARAMS_H_

#include "Solver/TimeStepping/Actor/ActorState.h"

#include <vector>

namespace seissol::solver {

/**
 * The part of a cluster's neighborhood that its step parameters depend on.
 *
 * Only configuration enters here (update rates and maximum time step sizes), never the progress
 * of a neighbor. That is what makes the step parameters independent of the order in which the
 * clusters are scheduled.
 */
struct StepContext {
  /// Update rate of the next-larger neighboring cluster, or 0 if there is none.
  long largerNeighborRate{0};

  /// Largest maximum time step size among the cluster and its neighbors.
  double largestTimeStepSize{0};
};

/**
 * Everything a single step of a cluster depends on, apart from the cell data.
 *
 * The step parameters are a function of the cluster's own times, the synchronization time and the
 * StepContext alone.
 *
 * On `subTimeStart`: a cluster with a next-larger neighbor integrates the derivatives of that
 * neighbor over its own step. These derivatives are expanded around the start of the neighbor's
 * step. Take a neighbor with a five times larger step:
 *
 *   neighbor:  |-------------------------------------------------------------|
 *   cluster:   |-----------|-----------|+++++++++++|-----------|-----------|
 *              0           dt          2dt         3dt         4dt         5dt
 *
 * For the step marked by `+`, the derivatives are integrated over [2dt, 3dt], i.e. `subTimeStart`
 * is 2dt. Since only the last step before a synchronization point is cut off, and that step also
 * ends the neighbor's step, all preceding steps inside the neighbor's step are complete; hence
 * `subTimeStart` is always a whole multiple of the cluster's maximum time step size.
 */
struct StepParams {
  /// Start time of the step.
  double time{0};

  /// Size of the step; the last step before the synchronization time is cut off there.
  double timeStepSize{0};

  /// Start of the step relative to the start of the enclosing step of the next-larger neighbor;
  /// 0 if there is no such neighbor.
  double subTimeStart{0};

  /// Time scale of the derivatives received from the next-larger neighbor.
  double neighborTimeStepSize{0};

  /// Whether the integrals accumulated for the next-larger neighbor start anew with this step.
  /// That is the case for the first step inside each step of the neighbor, and always if there is
  /// no such neighbor.
  bool resetBuffers{true};
};

/**
 * Extracts the StepContext of a cluster from its own times and its neighbors.
 */
StepContext computeStepContext(const ClusterTimes& times,
                               const std::vector<NeighborCluster>& neighbors);

/**
 * Computes the parameters of the step that starts at `times.correctionTime`.
 */
StepParams
    computeStepParams(const ClusterTimes& times, double syncTime, const StepContext& context);

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ACTOR_STEPPARAMS_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_FACECLUSTER_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_FACECLUSTER_H_

#include "Core/Common/Executor.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"

namespace seissol::solver {

/**
 * A cluster that works on faces between cells. The face work of a step needs the predictions of
 * the adjacent cells for that step, and the adjacent cells need its result for their correction.
 *
 * In terms of the actor model, the face work is the correction of the face cluster, while its
 * prediction does nothing. Since the face cluster declares its data ready only after its
 * correction, the adjacent cells wait for the face work before they correct.
 */
class FaceCluster : public ComputeCluster {
  public:
  ~FaceCluster() override = default;

  [[nodiscard]] DataReadiness dataReadiness() const final;

  /// Does the face work of one step on its own, outside of the time stepping (as
  /// CellCluster::predictOnce).
  void interactOnce(const StepParams& params) { interact(params); }

  protected:
  using ComputeCluster::ComputeCluster;

  /**
   * Does the face work of the step described by `params`.
   */
  virtual void interact(const StepParams& params) = 0;

  void start() override {}
  void predict() final {}
  void correct() final;

  /// the face work reads the cells (and the ghost layers) and writes the data of the faces; the
  /// prediction does nothing, so it neither waits nor has an event
  [[nodiscard]] DataAccess dataAccess(ActorAction action) const override {
    return action == ActorAction::Predict ? DataAccess::None : DataAccess::ReadWrite;
  }

  void handleNeighborPrediction(const NeighborCluster& /*neighborCluster*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighborCluster*/) override {}
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_FACECLUSTER_H_

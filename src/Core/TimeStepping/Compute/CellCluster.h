// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_

#include "Core/Common/Executor.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/CellStage.h"
#include "Core/TimeStepping/Compute/ClusterStage.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"

#include <functional>
#include <memory>
#include <vector>

namespace seissol::kernels {
class ReceiverCluster;
} // namespace seissol::kernels

namespace seissol::solver {

/**
 * The cells of one layer, i.e. of one time cluster and halo type, which all advance with the time
 * step of their time cluster: the frame of their steps. The numerics of the cells (see the clusters
 * of the wave propagation) fill in what a prediction and a correction compute; the frame runs that
 * work on the host or on the device, and the stages features add to the steps.
 **/
class CellCluster : public ComputeCluster {
  public:
  /**
   * The cells of `clusterData`, of the cluster `clusterId` on this rank (`globalClusterId` among
   * all), with the time step `maxTimeStepSize` and the rate `timeStepRate`, computing on
   * `executor`; `printProgress` makes it report its steps.
   **/
  CellCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              double maxTimeStepSize,
              long timeStepRate,
              Executor executor,
              bool printProgress,
              LTS::Layer* clusterData,
              LoopStatistics* loopStatistics,
              DeviceTimeline* deviceTimeline,
              ActorStateStatistics* actorStateStatistics);

  ~CellCluster() override = default;

  /// Adds work over items of its own, such as point sources, to the steps of the cluster.
  void addStage(std::unique_ptr<ClusterStage> stage);

  /// Adds work over the cells of the cluster to its corrections, such as the plasticity; the
  /// stages run in the order they were added.
  virtual void addCellStage(std::unique_ptr<CellStageBase> stage) = 0;

  /// Adds receivers the cluster samples; every output of receivers adds those of the cluster.
  virtual void addReceiverCluster(kernels::ReceiverCluster* receiverCluster) = 0;

  void finalize() override;

  [[nodiscard]] unsigned int getClusterId() const { return clusterId_; }
  [[nodiscard]] unsigned int getGlobalClusterId() const { return globalClusterId_; }

  protected:
  /// The work of a prediction of the cells on the host, as decided when it gets planned (e.g. the
  /// samples the receivers take); it runs later, on a thread of the host runtime.
  virtual std::function<void()> hostPrediction(const StepParams& params) = 0;
  /// The work of a prediction of the cells on the device, enqueued right away.
  virtual void devicePrediction(const StepParams& params) = 0;
  /// The work of a correction of the cells on the host, on a thread of the host runtime.
  virtual void hostCorrection(const StepParams& params) = 0;
  /// The work of a correction of the cells on the device, enqueued right away.
  virtual void deviceCorrection(const StepParams& params) = 0;

  /// Plans the prediction of the step `params` on the host, also for the device (e.g. counts its
  /// FLOPs); returns whether the step samples outputs. Only for a cluster with cells.
  virtual bool planPrediction(const StepParams& params) = 0;
  /// Plans the correction of the step `params` on the host, also for the device.
  virtual void planCorrection(const StepParams& params) = 0;

  /// Synchronizes what the neighbors of the cells read after a prediction (e.g. the buffers of the
  /// wave propagation) to `place`, on the stream of the cluster; for neighbors that compute on the
  /// other executor.
  virtual void synchronizeShared(initializer::AllocationPlace place) = 0;

  StepWork prepare(ActorAction action) override;
  void handleNeighborPrediction(const NeighborCluster& /*neighborCluster*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighborCluster*/) override {}
  void start() override {}
  void predict() override;
  void correct() override;

  LTS::Layer* clusterData_;

  private:
  /// Runs the stages of the phase `phase` for the step `params`.
  void runStages(StagePhase phase, const StepParams& params);

  //! the work features add to the steps, and the regions of the loop statistics they count in
  std::vector<std::unique_ptr<ClusterStage>> stages_;
  std::vector<unsigned> stageRegions_;

  //! report the steps
  bool printProgress_;
  //! cluster id on this rank
  unsigned int clusterId_;
  //! global cluster id
  unsigned int globalClusterId_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_

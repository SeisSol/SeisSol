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
#include "Core/Common/Real.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/ActorStateStatistics.h"
#include "Core/Monitoring/DeviceTimer.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/CellStage.h"
#include "Core/TimeStepping/Compute/ClusterClock.h"
#include "Core/TimeStepping/Compute/ClusterStage.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "WavePropagation/Output/FreeSurfaceIntegrator.h"
#include "WavePropagation/Output/Receiver.h"
#include "WavePropagation/Solver/ConfigBoundary.h"
#include "WavePropagation/Solver/Solver.h"
#include "WavePropagation/Solver/TimeCommon.h"

#include <cstdint>
#include <list>
#include <memory>
#include <mpi.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::kernels {
class ReceiverCluster;
} // namespace seissol::kernels

namespace seissol::solver {

/**
 * What the time manager sees of a cell cluster, whatever the configuration of its cells.
 **/
class CellClusterInterface : public ComputeCluster {
  public:
  using ComputeCluster::ComputeCluster;

  /// Adds work over items of its own, such as point sources, to the steps of the cluster.
  virtual void addStage(std::unique_ptr<ClusterStage> stage) = 0;

  /// Adds work over the cells of the cluster to its corrections, such as the plasticity; the
  /// stages run in the order they were added.
  virtual void addCellStage(std::unique_ptr<CellStageBase> stage) = 0;

  /// Adds receivers the cluster samples; every output of receivers adds those of the cluster.
  virtual void addReceiverCluster(kernels::ReceiverCluster* receiverCluster) = 0;
};

/**
 * The cells of one layer, i.e. of one time cluster and halo type, all of the configuration `Cfg`.
 * They all advance with the time step of their time cluster.
 **/
template <typename Cfg>
class CellCluster : public CellClusterInterface {
  public:
  using real = Real<Cfg>;

  private:
  seissol::SeisSol& seissolInstance_;
  /*
   * integrators
   */
  kernels::Spacetime<Cfg> spacetimeKernel_;
  //! time kernel
  kernels::Time<Cfg> timeKernel_;

  //! local kernel
  kernels::Local<Cfg> localKernel_;

  //! neighbor kernel
  kernels::Neighbor<Cfg> neighborKernel_;

  //! the neighbors of the cells that compute in other configurations
  kernels::ConfigBoundary<Cfg> configBoundary_;

  /*
   * global data
   */
  //! global data structures
  CompoundGlobalData<Cfg> globalData_;

  /*
   * element data
   */
  LTS::Layer* clusterData_;

  //! the work features add to the steps, and the regions of the loop statistics they count in
  std::vector<std::unique_ptr<ClusterStage>> stages_;
  std::vector<unsigned> stageRegions_;

  //! the work features add to the corrections over the cells, in the order it runs in; at most
  //! MaxCellStages of them, whose counts the loop over the cells sums up without allocating
  static constexpr std::size_t MaxCellStages = 4;
  std::vector<std::unique_ptr<CellStage<Cfg>>> cellStages_;

  enum class ComputePart : std::size_t { Local = 0, Neighbor, DRNeighbor, NumComputeParts };

  std::array<PerformanceEstimate, static_cast<std::size_t>(ComputePart::NumComputeParts)>
      estimate_{};
  std::array<std::size_t, static_cast<std::size_t>(ComputePart::NumComputeParts)> perfHandle_{};

  unsigned regionComputeLocalIntegration_;
  //! the names of the loops in a trace (see monitoring::Trace)
  std::uint32_t traceLocalIntegration_{};
  std::uint32_t traceNeighboringIntegration_{};
  unsigned regionComputeNeighboringIntegration_;

  //! the receivers of an output the cluster samples, the time of their next sample, and their
  //! samples in the current step, as prepare() plans them
  struct ReceiverSlot {
    kernels::ReceiverCluster* cluster{nullptr};
    double nextTime{0};
    kernels::ReceiverCluster::Sampling sampling{};
  };
  std::vector<ReceiverSlot> receivers_;
  //! the time the receivers start from (see setTime())
  double receiverStartTime_{0};

  /**
   * Samples the receivers of the cluster, with the samplings `samplings` planned for the step (one
   * per output, in the order of the outputs).
   **/
  void writeReceivers(const StepParams& params,
                      const std::vector<kernels::ReceiverCluster::Sampling>& samplings);

  /// the samplings of the receivers planned for the current step, in the order of the outputs
  [[nodiscard]] std::vector<kernels::ReceiverCluster::Sampling> plannedSamplings() const;

  /**
   * Runs the stages of the phase `phase` for the step `params`.
   **/
  void runStages(StagePhase phase, const StepParams& params);

  /**
   * Computes all cell local integration.
   *
   * This are:
   *  * time integration
   *  * volume integration
   *  * local boundary integration
   *
   * Remark: After this step the DOFs are only updated half with the boundary contribution
   *         of the neighborings cells missing.
   *
   * @param params parameters of the step.
   **/
  void computeLocalIntegration(const StepParams& params);

  /**
   * Computes the contribution of the neighboring cells to the boundary integral.
   *
   * Remark: After this step (in combination with the local integration) the DOFs are at the next
   *time step.
   * TODO: This excludes dynamic rupture contribution.
   *
   * @param params parameters of the step.
   **/
  void computeNeighboringIntegration(const StepParams& params);

  void computeLocalIntegrationDevice(const StepParams& params);
  void computeNeighboringIntegrationDevice(const StepParams& params);

  void computeLocalIntegrationFlops();

  PerformanceEstimate computeLocalIntegrationFlops(std::size_t numberOfCells,
                                                   const CellLocalInformation* cellInformation);

  void computeNeighborIntegrationFlops();

  void computeFlops();

  //! print status
  bool printProgress_;
  //! cluster id on this rank
  unsigned int clusterId_;

  //! global cluster cluster id
  unsigned int globalClusterId_;

  void incrementPerformanceMetrics(ComputePart part);

  protected:
  void handleNeighborPrediction(const NeighborCluster& neighborCluster) override;
  void handleNeighborCorrection(const NeighborCluster& neighborCluster) override;
  void start() override {}
  void predict() override;
  void correct() override;

  public:
  /**
   * Constructs a new LTS cluster.
   *
   * @param clusterId id of this cluster with respect to the current rank.
   * @param globalClusterId global id of this cluster.
   **/
  CellCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              double maxTimeStepSize,
              long timeStepRate,
              bool printProgress,
              CompoundGlobalData<Cfg> globalData,
              LTS::Layer* clusterData,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
              DeviceTimeline* deviceTimeline,
              ActorStateStatistics* actorStateStatistics);

  ~CellCluster() override = default;

  void addStage(std::unique_ptr<ClusterStage> stage) override;

  void addCellStage(std::unique_ptr<CellStageBase> stage) override;

  void addReceiverCluster(kernels::ReceiverCluster* receiverCluster) override {
    receiverCluster->setNextSampleTime(receiverStartTime_);
    receivers_.push_back(ReceiverSlot{.cluster = receiverCluster, .nextTime = receiverStartTime_});
  }

  void finalize() override;

  protected:
  StepWork prepare(ActorAction action) override;

  public:
  [[nodiscard]] bool outputsAhead(long steps) const override;
  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() override;

  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  void setTime(double time) override;

  void finishPhase() override;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_

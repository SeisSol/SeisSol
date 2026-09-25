// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_

#include "Common/Executor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Receiver.h"
#include "Kernels/Solver.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/Settings.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/StepParams.h"
#include "Solver/TimeStepping/Compute/ClusterClock.h"
#include "SourceTerm/Typedefs.h"

#include <list>
#include <memory>
#include <mpi.h>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::kernels {
class ReceiverCluster;
} // namespace seissol::kernels

namespace seissol::solver {

/**
 * The cells of one layer, i.e. of one time cluster and halo type. They all advance with the time
 * step of their time cluster.
 **/
class CellCluster : public AbstractTimeCluster {
  private:
  SimulationSettings settings_;

  seissol::SeisSol& seissolInstance_;
  /*
   * integrators
   */
  kernels::Spacetime spacetimeKernel_;
  //! time kernel
  kernels::Time timeKernel_;

  //! local kernel
  kernels::Local localKernel_;

  //! neighbor kernel
  kernels::Neighbor neighborKernel_;

  seissol::parallel::runtime::StreamRuntime streamRuntime_;

  //! the start of the current step, for the work on the stream
  ClusterClock clock_;

  /*
   * global data
   */
  //! global data structures
  GlobalData* globalDataOnHost_{nullptr};
  GlobalData* globalDataOnDevice_{nullptr};
#ifdef ACL_DEVICE
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif

  /*
   * element data
   */
  LTS::Layer* clusterData_;

  seissol::kernels::PointSourceClusterPair sourceCluster_;

  enum class ComputePart : std::size_t {
    Local = 0,
    Neighbor,
    DRNeighbor,
    PlasticityCheck,
    PlasticityYield,
    NumComputeParts
  };

  std::array<PerformanceEstimate, static_cast<std::size_t>(ComputePart::NumComputeParts)>
      estimate_{};
  std::array<std::size_t, static_cast<std::size_t>(ComputePart::NumComputeParts)> perfHandle_{};

  //! Stopwatch of TimeManager
  LoopStatistics* loopStatistics_;
  ActorStateStatistics* actorStateStatistics_;
  unsigned regionComputeLocalIntegration_;
  unsigned regionComputeNeighboringIntegration_;
  unsigned regionComputePointSources_;

  kernels::ReceiverCluster* receiverCluster_{nullptr};

  seissol::memory::MemkindArray<std::size_t> conditionalCounterHost_;
  seissol::memory::MemkindArray<std::size_t> conditionalCounterDevice_;

  std::size_t numPlasticCells_{0};

  /**
   * Writes the receiver output if applicable (receivers present, receivers have to be written).
   **/
  void writeReceivers(const StepParams& params);

  /**
   * Computes the source terms if applicable.
   **/
  void computeSources(const StepParams& params);

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

  template <bool UsePlasticity, bool IntegrateOutput>
  void computeNeighboringIntegrationImplementation(const StepParams& params);

  PerformanceEstimate computeLocalIntegrationFlops(std::size_t numberOfCells,
                                                   const CellLocalInformation* cellInformation);

  void computeNeighborIntegrationFlops();

  void computeFlops();

  HaloType layerType_;
  //! time of the next receiver output
  double receiverTime_;

  //! the receiver samples of the current step, as planned by prepare()
  kernels::ReceiverCluster::Sampling receiverSampling_{};

  //! the outputs decide about their samples when the work of a step runs
  bool runTimeOutputs_{false};

  //! print status
  bool printProgress_;
  //! cluster id on this rank
  unsigned int clusterId_;

  //! global cluster cluster id
  unsigned int globalClusterId_;

  //! id used to identify this cluster (including layer type) when profiling
  unsigned int profilingId_;

  void incrementPerformanceMetrics(ComputePart part);

  protected:
  void handleNeighborPrediction(const NeighborCluster& neighborCluster) override;
  void handleNeighborCorrection(const NeighborCluster& neighborCluster) override;
  void start() override {}
  void predict() override;
  void correct() override;

  public:
  ActResult act() override;

  /**
   * Constructs a new LTS cluster.
   *
   * @param clusterId id of this cluster with respect to the current rank.
   * @param globalClusterId global id of this cluster.
   * @param settings simulation settings
   **/
  CellCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              const SimulationSettings& settings,
              HaloType layerType,
              double maxTimeStepSize,
              long timeStepRate,
              bool printProgress,
              CompoundGlobalData globalData,
              LTS::Layer* clusterData,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
              ActorStateStatistics* actorStateStatistics);

  ~CellCluster() override = default;

  /**
   * Sets the the cluster's point sources
   *
   * @param sourceCluster Contains point sources for cluster
   */
  void setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster);

  void setReceiverCluster(kernels::ReceiverCluster* receiverCluster) {
    this->receiverCluster_ = receiverCluster;
    if (receiverCluster_ != nullptr) {
      receiverCluster_->setNextSampleTime(receiverTime_);
    }
  }

  void finalize() override;

  void setRunTimeOutputs(bool runTimeOutputs) override;

  protected:
  void* recordActionEvent() override;
  void waitForEvent(void* event) override;
  void timeSet(double time) override;
  StepWork prepare(ActorAction action) override;

  public:
  [[nodiscard]] bool outputsAhead(long steps) const override;
  [[nodiscard]] bool hostWork() const override;

  protected:
  public:
  [[nodiscard]] std::size_t layerId() const;
  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  [[nodiscard]] HaloType getLayerType() const;
  void setTime(double time) override;

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;

  void finishPhase() override;

  [[nodiscard]] std::string description() const override;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_CELLCLUSTER_H_

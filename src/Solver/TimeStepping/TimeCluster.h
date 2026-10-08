// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMECLUSTER_H_

#include "AbstractTimeCluster.h"
#include "Common/Executor.h"
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/ConfigBoundary.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Solver.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/Settings.h"
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

namespace seissol::time_stepping {

/**
 * What the time manager sees of a time cluster, whatever the configuration of its cells.
 **/
class TimeClusterInterface : public AbstractTimeCluster {
  public:
  using AbstractTimeCluster::AbstractTimeCluster;

  /**
   * Sets the the cluster's point sources
   *
   * @param sourceCluster Contains point sources for cluster
   */
  virtual void setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) = 0;

  virtual void setReceiverCluster(kernels::ReceiverCluster* receiverCluster) = 0;

  virtual void setFaultOutputManager(dr::output::OutputManager* outputManager) = 0;

  [[nodiscard]] virtual std::size_t layerId() const = 0;
};

/**
 * Time cluster, which represents a collection of elements having the same time step width, all of
 * the configuration `Cfg`.
 **/
template <typename Cfg>
class TimeCluster : public TimeClusterInterface {
  public:
  using real = Real<Cfg>; // NOLINT(readability-identifier-naming)

  private:
  // Last correction time of the neighboring cluster with higher dt
  double lastSubTime_{0};

  // The timestep of the largest neighbor. Not well-defined (and not used) for the largest local
  // timecluster.
  double neighborTimestep_{0};

  SimulationSettings settings_;

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

  kernels::DynamicRupture<Cfg> dynamicRuptureKernel_;

  seissol::parallel::runtime::StreamRuntime streamRuntime_;

  /*
   * global data
   */
  //! global data structures
  CompoundGlobalData<Cfg> globalData_;
#ifdef ACL_DEVICE
  device::DeviceInstance& device_ = device::DeviceInstance::instance();
#endif

  /*
   * element data
   */
  LTS::Layer* clusterData_;
  DynamicRupture::Layer* dynRupInteriorData_;
  DynamicRupture::Layer* dynRupCopyData_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolver_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolverDevice_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolverCopy_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolverCopyDevice_;
  dr::output::OutputManager* faultOutputManager_;

  seissol::kernels::PointSourceClusterPair sourceCluster_;

  enum class ComputePart : std::size_t {
    Local = 0,
    Neighbor,
    DRNeighbor,
    DRFrictionLawInterior,
    DRFrictionLawCopy,
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
  unsigned regionComputeDynamicRupture_;
  unsigned regionComputePointSources_;

  kernels::ReceiverCluster* receiverCluster_{nullptr};

  seissol::memory::MemkindArray<std::size_t> conditionalCounterHost_;
  seissol::memory::MemkindArray<std::size_t> conditionalCounterDevice_;

  std::size_t numPlasticCells_{0};

  /**
   * Writes the receiver output if applicable (receivers present, receivers have to be written).
   **/
  void writeReceivers();

  /**
   * Computes the source terms if applicable.
   **/
  void computeSources();

  /**
   * Computes dynamic rupture.
   **/
  void computeDynamicRupture(DynamicRupture::Layer& layerData);

  void handleDynamicRupture(DynamicRupture::Layer& layerData);

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
   * @param numberOfCells number of cells.
   * @param cellInformation cell local information.
   * @param cellData cell data.
   * @param buffers time integration buffers.
   * @param derivatives time derivatives.
   * @param dofs degrees of freedom.
   **/
  void computeLocalIntegration(bool resetBuffers);

  /**
   * Computes the contribution of the neighboring cells to the boundary integral.
   *
   * Remark: After this step (in combination with the local integration) the DOFs are at the next
   *time step.
   * TODO: This excludes dynamic rupture contribution.
   *
   * @param numberOfCells number of cells.
   * @param cellInformation cell local information.
   * @param cellData cell data.
   * @param faceNeighbors pointers to neighboring time buffers or derivatives.
   * @param dofs degrees of freedom.
   **/
  void computeNeighboringIntegration(double subTimeStart);

  void computeLocalIntegrationDevice(bool resetBuffers);
  void computeDynamicRuptureDevice(DynamicRupture::Layer& layerData);
  void computeNeighboringIntegrationDevice(double subTimeStart);

  void computeLocalIntegrationFlops();

  template <bool UsePlasticity, bool IntegrateOutput>
  void computeNeighboringIntegrationImplementation(double subTimeStart);

  PerformanceEstimate computeLocalIntegrationFlops(std::size_t numberOfCells,
                                                   const CellLocalInformation* cellInformation);

  void computeNeighborIntegrationFlops();

  PerformanceEstimate computeDynamicRuptureFlops(DynamicRupture::Layer& layerData);

  void computeFlops();

  HaloType layerType_;
  //! time of the next receiver output
  double receiverTime_;

  //! print status
  bool printProgress_;
  //! cluster id on this rank
  unsigned int clusterId_;

  //! global cluster cluster id
  unsigned int globalClusterId_;

  //! id used to identify this cluster (including layer type) when profiling
  unsigned int profilingId_;

  DynamicRuptureScheduler* dynamicRuptureScheduler_;

  void incrementPerformanceMetrics(ComputePart part);

  protected:
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
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
  TimeCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              const SimulationSettings& settings,
              HaloType layerType,
              double maxTimeStepSize,
              long timeStepRate,
              bool printProgress,
              DynamicRuptureScheduler* dynamicRuptureScheduler,
              CompoundGlobalData<Cfg> globalData,
              LTS::Layer* clusterData,
              DynamicRupture::Layer* dynRupInteriorData,
              DynamicRupture::Layer* dynRupCopyData,
              const seissol::dr::friction_law::FrictionSolverFactory& frictionSolverFactory,
              const seissol::dr::friction_law::FrictionSolverFactory& frictionSolverFactoryDevice,
              dr::output::OutputManager* faultOutputManager,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
              ActorStateStatistics* actorStateStatistics);

  ~TimeCluster() override = default;

  void setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) override;

  void setReceiverCluster(kernels::ReceiverCluster* receiverCluster) override {
    this->receiverCluster_ = receiverCluster;
  }

  void setFaultOutputManager(dr::output::OutputManager* outputManager) override {
    faultOutputManager_ = outputManager;
  }

  void reset() override;

  void finalize() override;

  [[nodiscard]] std::size_t layerId() const override;
  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  [[nodiscard]] HaloType getLayerType() const;
  void setTime(double time) override;

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;

  void finishPhase() override;

  [[nodiscard]] std::string description() const override;

  /**
   * A time step that equals the cluster's own step size recurs for the whole run, so a compute
   * graph recorded for it pays off. The truncated step taken right before a synchronization
   * point does not recur: caching it would add one graph per synchronization point.
   */
  [[nodiscard]] bool isRecurringTimestep(double timestep) const {
    return timestep == ct_.maxTimeStepSize;
  }
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMECLUSTER_H_

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

#include <list>
#include <mpi.h>

#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "SourceTerm/Typedefs.h"
#include <utils/logger.h>

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Solver.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/LoopStatistics.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include <Common/Executor.h>

#include "AbstractTimeCluster.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::kernels {
class ReceiverCluster;
} // namespace seissol::kernels

namespace seissol::time_stepping {

/**
 * Time cluster, which represents a collection of elements having the same time step width.
 **/
class TimeCluster : public AbstractTimeCluster {
  private:
  // Last correction time of the neighboring cluster with higher dt
  double lastSubTime;

  // The timestep of the largest neighbor. Not well-defined (and not used) for the largest local
  // timecluster.
  double neighborTimestep;

  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void start() override {}
  void predict() override;
  void correct() override;
  bool usePlasticity;

  //! number of time steps
  unsigned long m_numberOfTimeSteps;

  seissol::SeisSol& seissolInstance;
  /*
   * integrators
   */
  kernels::Spacetime spacetimeKernel;
  //! time kernel
  kernels::Time timeKernel;

  //! local kernel
  kernels::Local localKernel;

  //! neighbor kernel
  kernels::Neighbor neighborKernel;

  kernels::DynamicRupture dynamicRuptureKernel;

  seissol::parallel::runtime::StreamRuntime streamRuntime;

  /*
   * global data
   */
  //! global data structures
  GlobalData* globalDataOnHost{nullptr};
  GlobalData* globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  /*
   * element data
   */
  LTS::Layer* clusterData;
  DynamicRupture::Layer* dynRupInteriorData;
  DynamicRupture::Layer* dynRupCopyData;
  dr::friction_law::FrictionSolver* frictionSolver;
  dr::friction_law::FrictionSolver* frictionSolverDevice;
  dr::output::OutputManager* faultOutputManager;

  seissol::kernels::PointSourceClusterPair sourceCluster;

  enum class ComputePart : int {
    Local = 0,
    Neighbor,
    DRNeighbor,
    DRFrictionLawInterior,
    DRFrictionLawCopy,
    PlasticityCheck,
    PlasticityYield,
    NumComputeParts
  };

  std::array<std::uint64_t, static_cast<int>(ComputePart::NumComputeParts)> accFlopsNonZero;
  std::array<std::uint64_t, static_cast<int>(ComputePart::NumComputeParts)> accFlopsHardware;

  //! Stopwatch of TimeManager
  LoopStatistics* loopStatistics;
  ActorStateStatistics* actorStateStatistics;
  unsigned regionComputeLocalIntegration;
  unsigned regionComputeNeighboringIntegration;
  unsigned regionComputeDynamicRupture;
  unsigned regionComputePointSources;

  kernels::ReceiverCluster* receiverCluster{nullptr};

  seissol::memory::MemkindArray<std::size_t> yieldCells;

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

#ifdef ACL_DEVICE
  void computeLocalIntegrationDevice(bool resetBuffers);
  void computeDynamicRuptureDevice(DynamicRupture::Layer& layerData);
  void computeNeighboringIntegrationDevice(double subTimeStart);
#endif

  void computeLocalIntegrationFlops();

  template <bool UsePlasticity>
  void computeNeighboringIntegrationImplementation(double subTimeStart);

  void computeLocalIntegrationFlops(unsigned numberOfCells,
                                    const CellLocalInformation* cellInformation,
                                    std::uint64_t& nonZeroFlops,
                                    std::uint64_t& hardwareFlops);

  void computeNeighborIntegrationFlops();

  void computeDynamicRuptureFlops(DynamicRupture::Layer& layerData,
                                  std::uint64_t& nonZeroFlops,
                                  std::uint64_t& hardwareFlops);

  void computeFlops();

  const LayerType layerType;
  //! time of the next receiver output
  double receiverTime;

  //! print status every 100th timestep
  bool printProgress;
  //! cluster id on this rank
  const unsigned int clusterId;

  //! global cluster cluster id
  const unsigned int globalClusterId;

  //! id used to identify this cluster (including layer type) when profiling
  const unsigned int profilingId;

  DynamicRuptureScheduler* dynamicRuptureScheduler;

  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  public:
  ActResult act() override;

  /**
   * Constructs a new LTS cluster.
   *
   * @param clusterId id of this cluster with respect to the current rank.
   * @param globalClusterId global id of this cluster.
   * @param usePlasticity true if using plasticity
   **/
  TimeCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              bool usePlasticity,
              LayerType layerType,
              double maxTimeStepSize,
              long timeStepRate,
              bool printProgress,
              DynamicRuptureScheduler* dynamicRuptureScheduler,
              CompoundGlobalData globalData,
              LTS::Layer* clusterData,
              DynamicRupture::Layer* dynRupInteriorData,
              DynamicRupture::Layer* dynRupCopyData,
              seissol::dr::friction_law::FrictionSolver* frictionSolver,
              seissol::dr::friction_law::FrictionSolver* frictionSolverDevice,
              dr::output::OutputManager* faultOutputManager,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
              ActorStateStatistics* actorStateStatistics);

  ~TimeCluster() override = default;

  /**
   * Sets the the cluster's point sources
   *
   * @param sourceCluster Contains point sources for cluster
   */
  void setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster);

  void setReceiverCluster(kernels::ReceiverCluster* receiverCluster) {
    this->receiverCluster = receiverCluster;
  }

  void setFaultOutputManager(dr::output::OutputManager* outputManager) {
    faultOutputManager = outputManager;
  }

  void reset() override;

  void finalize() override;

  [[nodiscard]] std::size_t layerId() const;
  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  [[nodiscard]] LayerType getLayerType() const;
  void setTime(double time) override;

  std::vector<NeighborCluster>* getNeighborClusters();

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;

  void finishPhase() override;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMECLUSTER_H_

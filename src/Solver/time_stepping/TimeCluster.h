// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_TIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_TIMECLUSTER_H_

#ifdef USE_MPI
#include <mpi.h>
#include <list>
#endif

#include "Initializer/Typedefs.h"
#include "SourceTerm/Typedefs.h"
#include <utils/logger.h>
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"

#include "Kernels/Time.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/TimeCommon.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include <Common/Executor.h>

#include "AbstractTimeCluster.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
  namespace time_stepping {
    class TimeCluster;
  }

  namespace kernels {
    class ReceiverCluster;
  }
}

/**
 * Time cluster, which represents a collection of elements having the same time step width.
 **/
class seissol::time_stepping::TimeCluster : public seissol::time_stepping::AbstractTimeCluster
{
private:
    // Last correction time of the neighboring cluster with higher dt
    double lastSubTime;

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
    //! time kernel
    kernels::Time m_timeKernel;

    //! local kernel
    kernels::Local m_localKernel;

    //! neighbor kernel
    kernels::Neighbor m_neighborKernel;
    
    kernels::DynamicRupture m_dynamicRuptureKernel;

    seissol::parallel::runtime::StreamRuntime streamRuntime;

  /*
   * global data
   */
     //! global data structures
    GlobalData *m_globalDataOnHost{nullptr};
    GlobalData *m_globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

    /*
     * element data
     */     
    seissol::initializer::Layer* m_clusterData;
    seissol::initializer::Layer* dynRupInteriorData;
    seissol::initializer::Layer* dynRupCopyData;
    seissol::initializer::LTS*         m_lts;
    seissol::initializer::DynamicRupture* m_dynRup;
    dr::friction_law::FrictionSolver* frictionSolver;
    dr::friction_law::FrictionSolver* frictionSolverDevice;
    dr::output::OutputManager* faultOutputManager;

    seissol::kernels::PointSourceClusterPair m_sourceCluster;

    enum class ComputePart {
      Local = 0,
      Neighbor,
      DRNeighbor,
      DRFrictionLawInterior,
      DRFrictionLawCopy,
      PlasticityCheck,
      PlasticityYield,
      NUM_COMPUTE_PARTS
    };

    long long m_flops_nonZero[static_cast<int>(ComputePart::NUM_COMPUTE_PARTS)];
    long long m_flops_hardware[static_cast<int>(ComputePart::NUM_COMPUTE_PARTS)];
    
    //! Tv parameter for plasticity
    double m_tv;
    
    //! Relax time for plasticity
    double m_oneMinusIntegratingFactor;
    
    //! Stopwatch of TimeManager
    LoopStatistics* m_loopStatistics;
    ActorStateStatistics* actorStateStatistics;
    unsigned        m_regionComputeLocalIntegration;
    unsigned        m_regionComputeNeighboringIntegration;
    unsigned        m_regionComputeDynamicRupture;
    unsigned        m_regionComputePointSources;

    kernels::ReceiverCluster* m_receiverCluster;

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
    void computeDynamicRupture( seissol::initializer::Layer&  layerData );

    void handleDynamicRupture( seissol::initializer::Layer&  layerData );

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
     * @param i_numberOfCells number of cells.
     * @param i_cellInformation cell local information.
     * @param i_cellData cell data.
     * @param io_buffers time integration buffers.
     * @param io_derivatives time derivatives.
     * @param io_dofs degrees of freedom.
     **/
    void computeLocalIntegration( seissol::initializer::Layer&  layerData, bool resetBuffers);

    /**
     * Computes the contribution of the neighboring cells to the boundary integral.
     *
     * Remark: After this step (in combination with the local integration) the DOFs are at the next time step.
     * TODO: This excludes dynamic rupture contribution.
     *
     * @param i_numberOfCells number of cells.
     * @param i_cellInformation cell local information.
     * @param i_cellData cell data.
     * @param i_faceNeighbors pointers to neighboring time buffers or derivatives.
     * @param io_dofs degrees of freedom.
     **/
    void computeNeighboringIntegration( seissol::initializer::Layer&  layerData, double subTimeStart );

#ifdef ACL_DEVICE
    void computeLocalIntegrationDevice( seissol::initializer::Layer&  layerData, bool resetBuffers);
    void computeDynamicRuptureDevice( seissol::initializer::Layer&  layerData );
    void computeNeighboringIntegrationDevice( seissol::initializer::Layer&  layerData, double subTimeStart );
#endif

    void computeLocalIntegrationFlops(seissol::initializer::Layer& layerData);

    template<bool usePlasticity>
    std::pair<long, long> computeNeighboringIntegrationImplementation(seissol::initializer::Layer& layerData,
                                                                      double subTimeStart);

    void computeLocalIntegrationFlops(unsigned numberOfCells,
                                      CellLocalInformation const* cellInformation,
                                      long long& nonZeroFlops,
                                      long long& hardwareFlops);

    void computeNeighborIntegrationFlops(seissol::initializer::Layer &layerData);

    void computeDynamicRuptureFlops(seissol::initializer::Layer &layerData,
                                    long long& nonZeroFlops,
                                    long long& hardwareFlops);
                                          
    void computeFlops();
    
    //! Update relax time for plasticity
    void updateRelaxTime() {
      m_oneMinusIntegratingFactor = (m_tv > 0.0) ? 1.0 - exp(-timeStepSize() / m_tv) : 1.0;
    }

  const LayerType layerType;
  //! time of the next receiver output
  double m_receiverTime;

  //! print status every 100th timestep
  bool printProgress;
  //! cluster id on this rank
  const unsigned int m_clusterId;

  //! global cluster cluster id
  const unsigned int m_globalClusterId;

  //! id used to identify this cluster (including layer type) when profiling
  const unsigned int m_profilingId;

  DynamicRuptureScheduler* dynamicRuptureScheduler;

  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

public:
  ActResult act() override;

  /**
   * Constructs a new LTS cluster.
   *
   * @param i_clusterId id of this cluster with respect to the current rank.
   * @param i_globalClusterId global id of this cluster.
   * @param usePlasticity true if using plasticity
   **/
  TimeCluster(unsigned int i_clusterId,
      unsigned int i_globalClusterId,
      unsigned int profilingId,
      bool usePlasticity,
      LayerType layerType,
      double maxTimeStepSize,
      long timeStepRate,
      bool printProgress,
      DynamicRuptureScheduler* dynamicRuptureScheduler,
      CompoundGlobalData i_globalData,
      seissol::initializer::Layer *i_clusterData,
      seissol::initializer::Layer* dynRupInteriorData,
      seissol::initializer::Layer* dynRupCopyData,
      seissol::initializer::LTS* i_lts,
      seissol::initializer::DynamicRupture* i_dynRup,
      seissol::dr::friction_law::FrictionSolver* i_FrictionSolver,
      seissol::dr::friction_law::FrictionSolver* i_FrictionSolverDevice,
      dr::output::OutputManager* i_faultOutputManager,
      seissol::SeisSol& seissolInstance,
      LoopStatistics* i_loopStatistics,
      ActorStateStatistics* actorStateStatistics);

  /**
   * Destructor of a LTS cluster.
   * TODO: Currently prints only statistics in debug mode.
   **/
  ~TimeCluster() override;

  /**
   * Sets the the cluster's point sources
   *
   * @param sourceCluster Contains point sources for cluster
   */
  void setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster);
  void freePointSources() { m_sourceCluster.host.reset(nullptr); m_sourceCluster.device.reset(nullptr); }

  void setReceiverCluster( kernels::ReceiverCluster* receiverCluster) {
    m_receiverCluster = receiverCluster;
  }

  void setFaultOutputManager(dr::output::OutputManager* outputManager) {
    faultOutputManager = outputManager;
  }

  /**
   * Set Tv constant for plasticity.
   */
  void setTv(double tv) {
    m_tv = tv;
    updateRelaxTime();
  }

  void setLastSubTime(double lastSubTime) {
    this->lastSubTime = lastSubTime;
  }


  void reset() override;

  void finalize() override;

  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  [[nodiscard]] LayerType getLayerType() const;
  void setReceiverTime(double receiverTime);

  std::vector<NeighborCluster>* getNeighborClusters();

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream);
};


#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_TIMECLUSTER_H_


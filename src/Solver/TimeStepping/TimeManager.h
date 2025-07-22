// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMEMANAGER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMEMANAGER_H_
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <cassert>
#include <list>
#include <memory>
#include <queue>
#include <vector>

#include "Initializer/MemoryManager.h"
#include "Initializer/TimeStepping/LtsLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/PointSourceCluster.h"
#include "Monitoring/Stopwatch.h"
#include "ResultWriter/ReceiverWriter.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/TimeStepping/GhostTimeClusterFactory.h"
#include "SourceTerm/Typedefs.h"
#include "TimeCluster.h"
#include <utils/logger.h>

namespace seissol::time_stepping {
class AbstractCommunicationManager;

/**
 * Time manager, which takes care of the time stepping.
 **/
class TimeManager {
  private:
  seissol::SeisSol& seissolInstance;

  //! time stepping
  std::optional<initializer::ClusterLayout> clusterLayout;

  //! all local (copy & interior) LTS clusters, which are under control of this time manager
  std::vector<std::unique_ptr<TimeCluster>> clusters;
  std::vector<TimeCluster*> highPrioClusters;
  std::vector<TimeCluster*> lowPrioClusters;

  //! one dynamic rupture scheduler per pair of interior/copy cluster
  std::vector<std::unique_ptr<DynamicRuptureScheduler>> dynamicRuptureSchedulers;

  //! all MPI (ghost) LTS clusters, which are under control of this time manager
  std::unique_ptr<AbstractCommunicationManager> communicationManager;

  //! Stopwatch
  LoopStatistics loopStatistics;
  ActorStateStatisticsManager actorStateStatisticsManager;

  //! dynamic rupture output
  dr::output::OutputManager* faultOutputManager{};

  public:
  /**
   * Construct a new time manager.
   **/
  TimeManager(seissol::SeisSol& seissolInstance);

  /**
   * Destruct the time manager.
   **/
  ~TimeManager();

  /**
   * Adds the time clusters to the time manager.
   *
   * @param i_timeStepping time stepping scheme.
   * @param i_meshStructure mesh structure.
   * @param memoryManager memory manager.
   * @param i_meshToClusters mapping from the mesh to the clusters.
   **/
  void addClusters(const initializer::ClusterLayout& clusterLayout,
                   const solver::HaloCommunication& haloStructure,
                   initializer::MemoryManager& memoryManager,
                   bool usePlasticity);

  void setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager);
  seissol::dr::output::OutputManager* getFaultOutputManager();

  /**
   * Advance in time until all clusters reach the next synchronization time.
   **/
  void advanceInTime(const double& synchronizationTime);

  /**
   * Gets the time tolerance of the time manager (1E-5 of the CFL time step width).
   **/
  [[nodiscard]] double getTimeTolerance() const;

  /**
   * Distributes point sources pointers to clusters
   *
   * @param sourceClusters Collection of point sources for clusters
   */
  void setPointSourcesForClusters(
      std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>>
          sourceClusters);

  /**
   * Returns the writer for the receivers
   */
  void setReceiverClusters(writer::ReceiverWriter& receiverWriter);

  /**
   * Sets the initial time (time DOFS/DOFs/receivers) of all time clusters.
   * Required only if different from zero, for example in checkpointing.
   *
   * @param i_time time.
   **/
  void setInitialTimes(double time = 0);

  void printComputationTime(const std::string& outputPrefix, bool isLoopStatisticsNetcdfOutputOn);

  void freeDynamicResources();

  void synchronizeTo(seissol::initializer::AllocationPlace place);

  const initializer::ClusterLayout& getClusterLayout() { return clusterLayout.value(); }
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMEMANAGER_H_

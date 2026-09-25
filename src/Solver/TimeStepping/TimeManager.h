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
#include "CellCluster.h"
#include "DynamicRuptureCluster.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Kernels/PointSourceCluster.h"
#include "Monitoring/Stopwatch.h"
#include "ResultWriter/ReceiverWriter.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "Solver/TimeStepping/GhostCluster.h"
#include "Solver/TimeStepping/HaloTransport.h"
#include "Solver/TimeStepping/SuperStepRecorder.h"
#include "Solver/TimeStepping/TimeSteppingPlan.h"
#include "SourceTerm/Typedefs.h"

#include <cassert>
#include <list>
#include <memory>
#include <queue>
#include <set>
#include <utils/logger.h>
#include <vector>

namespace seissol::time_stepping {
class AbstractCommunicationManager;

/**
 * Time manager, which takes care of the time stepping.
 **/
class TimeManager {
  private:
  seissol::SeisSol& seissolInstance_;

  //! time stepping
  std::optional<initializer::ClusterLayout> clusterLayout_;

  //! the clusters of the local (copy & interior) cells
  std::vector<std::unique_ptr<CellCluster>> cellClusters_;

  //! the clusters of the local (copy & interior) dynamic rupture faces
  std::vector<std::unique_ptr<DynamicRuptureCluster>> faceClusters_;

  //! all cell and face clusters, i.e. all clusters under control of this time manager that do
  //! local work; ordered by rate
  std::vector<AbstractTimeCluster*> clusters_;
  std::vector<AbstractTimeCluster*> highPrioClusters_;
  std::vector<AbstractTimeCluster*> lowPrioClusters_;

  //! what the halo transports share; outlives them
  std::unique_ptr<HaloTransportFactory> haloTransports_;

  //! all MPI (ghost) LTS clusters, which are under control of this time manager
  std::unique_ptr<AbstractCommunicationManager> communicationManager_;

  //! take the steps along the time stepping plan
  bool followPlan_{false};

  //! along the plan: all super-timesteps (steps of the largest cluster), those that end early at a
  //! synchronization point, the full ones without output samples, and the full ones in which no
  //! cluster takes an irregular step
  std::size_t superSteps_{0};
  std::size_t shortenedSuperSteps_{0};
  std::size_t outputFreeSuperSteps_{0};
  std::size_t regularSuperSteps_{0};
  std::size_t recordedSuperSteps_{0};
  std::size_t replayedSuperSteps_{0};
  std::size_t mispredictedSuperSteps_{0};

  //! record regular super-timesteps into graphs and replay them
  bool replay_{false};
  std::unique_ptr<SuperStepRecorder> recorder_;
  //! the regular super-timesteps that have run once without a recording
  std::set<SuperStepRecorder::Key> seenSuperSteps_;

  /**
   * Takes the steps [begin, end) of the plan, and returns what their host parts have decided.
   */
  StepWork takeSteps(const std::vector<PlannedAction>& plan, std::size_t begin, std::size_t end);

  //! the clusters only enqueue their device work
  bool concurrent_{false};

  /**
   * Takes all steps of the cell and face clusters up to the synchronization point in the order of
   * the time stepping plan.
   */
  void followPlan();

  //! Stopwatch
  LoopStatistics loopStatistics_;
  ActorStateStatisticsManager actorStateStatisticsManager_;

  //! dynamic rupture output
  dr::output::OutputManager* faultOutputManager_{};

  public:
  /**
   * Construct a new time manager.
   **/
  explicit TimeManager(seissol::SeisSol& seissolInstance);

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
                   const SimulationSettings& settings);

  void setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager);
  seissol::dr::output::OutputManager* faultOutputManager();

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
      std::vector<seissol::kernels::PointSourceClusterPair> sourceClusters);

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

  const initializer::ClusterLayout& getClusterLayout() { return clusterLayout_.value(); }
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_TIMEMANAGER_H_

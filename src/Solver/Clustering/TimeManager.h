// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 **/

#ifndef TIMEMANAGER_H_
#define TIMEMANAGER_H_
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Parallel/Host/CpuExecutor.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
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
#include "Solver/Clustering/Communication/CommunicationManager.h"
#include "Solver/Clustering/Communication/GhostTimeClusterFactory.h"
#include "Solver/Clustering/Computation/DynamicRuptureCluster.h"
#include "Solver/Clustering/Computation/TimeCluster.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "SourceTerm/Typedefs.h"
#include <utils/logger.h>

namespace seissol::solver::clustering {

template <typename T>
constexpr T ipow(T x, T y) {
  static_assert(std::is_integral_v<T>);
  assert(y >= 0);

  if (y == 0) {
    return 1;
  }
  T result = x;
  while (--y) {
    result *= x;
  }
  return result;
}

/**
 * Time manager, which takes care of the time stepping.
 **/
class TimeManager {
  private:
  //! last #updates of log
  unsigned int logUpdates;

  seissol::SeisSol& seissolInstance;

  //! time stepping
  TimeStepping timeStepping;

  std::shared_ptr<parallel::host::CpuExecutor> cpuExecutor;

  //! all local (copy & interior) LTS clusters, which are under control of this time manager
  // std::vector<std::unique_ptr<AbstractTimeCluster>> clusters;
  // std::vector<std::weak_ptr<AbstractTimeCluster>> highPrioClusters;
  // std::vector<std::weak_ptr<AbstractTimeCluster>> lowPrioClusters;

  std::vector<std::shared_ptr<AbstractTimeCluster>> clusters;

  std::vector<computation::TimeCluster*> highPrioClusters;
  std::vector<computation::TimeCluster*> lowPrioClusters;
  std::vector<computation::DynamicRuptureCluster*> highPrioClustersDR;
  std::vector<computation::DynamicRuptureCluster*> lowPrioClustersDR;

  //! all MPI (ghost) LTS clusters, which are under control of this time manager
  std::unique_ptr<communication::AbstractCommunicationManager> communicationManager;

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

  void addClusters(initializer::ClusterLayout& layout,
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
  double getTimeTolerance();

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
   * Set Tv constant for plasticity.
   */
  void setTv(double tv);

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

  inline const TimeStepping* getTimeStepping() { return &timeStepping; }
};

} // namespace seissol::solver::clustering

#endif

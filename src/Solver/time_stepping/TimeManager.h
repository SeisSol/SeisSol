// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_TIMEMANAGER_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_TIMEMANAGER_H_
#include <vector>
#include <queue>
#include <list>
#include <cassert>
#include <memory>

#include "Initializer/Typedefs.h"
#include "SourceTerm/Typedefs.h"
#include <utils/logger.h>
#include "Initializer/MemoryManager.h"
#include "Initializer/TimeStepping/LtsLayout.h"
#include "Kernels/PointSourceCluster.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include "ResultWriter/ReceiverWriter.h"
#include "TimeCluster.h"
#include "Monitoring/Stopwatch.h"
#include "Solver/time_stepping/GhostTimeClusterFactory.h"

namespace seissol {
  namespace time_stepping {
    class TimeManager;
    class AbstractCommunicationManager;

      template<typename T>
      constexpr T ipow(T x, T y) {
          static_assert(std::is_integral_v<T>);
          assert(y >= 0);

          if (y == 0) {
              return 1;
          }
          T result = x;
          while(--y) {
              result *= x;
          }
          return result;
      }
  }
}


/**
 * Time manager, which takes care of the time stepping.
 **/
class seissol::time_stepping::TimeManager {
  //private:
    /**
     * Compares to cluster pointer by their id.
     **/
    struct clusterCompare {
      bool operator()( const TimeCluster* l_first, const TimeCluster* l_second ) {
        return l_first->getGlobalClusterId() > l_second->getGlobalClusterId();
      }
    };

    //! last #updates of log
    unsigned int m_logUpdates;

    seissol::SeisSol& seissolInstance;

    //! time stepping
    TimeStepping m_timeStepping;

    //! all local (copy & interior) LTS clusters, which are under control of this time manager
    std::vector<std::unique_ptr<TimeCluster>> clusters;
    std::vector<TimeCluster*> highPrioClusters;
    std::vector<TimeCluster*> lowPrioClusters;

    //! one dynamic rupture scheduler per pair of interior/copy cluster
    std::vector<std::unique_ptr<DynamicRuptureScheduler>> dynamicRuptureSchedulers;

    //! all MPI (ghost) LTS clusters, which are under control of this time manager
    std::unique_ptr<AbstractCommunicationManager> communicationManager;

    //! Stopwatch
    LoopStatistics m_loopStatistics;
    ActorStateStatisticsManager actorStateStatisticsManager;
    
    //! dynamic rupture output
    dr::output::OutputManager* m_faultOutputManager{};

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
    void addClusters(TimeStepping& i_timeStepping,
                     MeshStructure* i_meshStructure,
                     initializer::MemoryManager& memoryManager,
                     bool usePlasticity);

    void setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager);
    seissol::dr::output::OutputManager* getFaultOutputManager();

    /**
     * Advance in time until all clusters reach the next synchronization time.
     **/
    void advanceInTime( const double &synchronizationTime );

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
        std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>> sourceClusters);

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
    void setInitialTimes( double i_time = 0 );

    void printComputationTime(const std::string& outputPrefix, bool isLoopStatisticsNetcdfOutputOn);

    void freeDynamicResources();

    void synchronizeTo(seissol::initializer::AllocationPlace place);

    inline const TimeStepping* getTimeStepping() {
      return &m_timeStepping;
    }
};


#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_TIMEMANAGER_H_


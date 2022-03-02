/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 * 
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Time step width management in SeisSol.
 **/

#ifndef TIMEMANAGER_H_
#define TIMEMANAGER_H_
#include <vector>
#include <queue>
#include <list>
#include <cassert>
#include <memory>

#include <Initializer/typedefs.hpp>
#include <SourceTerm/typedefs.hpp>
#include <utils/logger.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/time_stepping/LtsLayout.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <ResultWriter/ReceiverWriter.h>
#include "TimeCluster.h"
#include "Monitoring/Stopwatch.h"
#include "GhostTimeCluster.h"

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
    
  public:
    /**
     * Construct a new time manager.
     **/
    TimeManager();

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
                     initializers::MemoryManager& memoryManager,
                     bool usePlasticity);

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
     * @param clusterMappings Maps layers+clusters to point sources
     * @param pointSources Map from layer to list of point sources
     */
    void setPointSourcesForClusters(
        std::unordered_map<LayerType, std::vector<sourceterm::ClusterMapping>>& clusterMappings,
        std::unordered_map<LayerType, std::vector<sourceterm::PointSources>>& pointSources);

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

    void printComputationTime();
};

#endif

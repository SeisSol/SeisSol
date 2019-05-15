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

#include <Initializer/typedefs.hpp>
#include <SourceTerm/typedefs.hpp>
#include <utils/logger.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/time_stepping/LtsLayout.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <ResultWriter/ReceiverWriter.h>
#include "TimeCluster.h"
#include "Monitoring/Stopwatch.h"

namespace seissol {
  namespace time_stepping {
    class TimeManager;
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
        return l_first->m_globalClusterId > l_second->m_globalClusterId;
      }
    };

    //! last #updates of log
    unsigned int m_logUpdates;

    //! time stepping
    TimeStepping m_timeStepping;

    //! all LTS clusters, which are under control of this time manager
    std::vector< TimeCluster* > m_clusters;

    //! queue of clusters, which are allowed to update their copy layer locally
    std::list< TimeCluster* > m_localCopyQueue;

    //! queue of clusters, which are allowed to update their interior locally
    std::priority_queue< TimeCluster*, std::vector<TimeCluster*>, clusterCompare > m_localInteriorQueue;

    //! queue of clusters which are allowed to update their copy layer with the neighboring cells contribution
    std::list< TimeCluster* > m_neighboringCopyQueue;

    //! queue of clusters which are allowed to update their interior with the neighboring cells contribution
    std::priority_queue< TimeCluster*, std::vector<TimeCluster*>, clusterCompare > m_neighboringInteriorQueue;
    
    //! Stopwatch
    LoopStatistics m_loopStatistics;
    
    /**
     * Checks if the time stepping restrictions for this cluster and its neighbors changed.
     * If this is true:
     *  * The respectice queues are updated.
     *  * The corresponding copy layer/interior is set eligible for an full update or prediction.
     *  * In the case of a new prediction the next time step width is derived.
     *  * Enables the reset of the time buffers if applicable.
     *
     * Sketch:
     *   ________ ______ ________
     *  |        |      |        |
     *  | TC n-1 | TC n | TC n+1 |
     *  |________|______|________|
     *
     * A time cluster with id n is a neighbor of time clusters n-1 and n+1.
     * By performing a full time step update or providing a new time predicion in cluster n
     * clusters n-1 or n or n+1 migh now be allowed to perform a full update or a new prediction.
     *
     * @param i_localClusterId local cluster id of the cluster, which changed its status.
     **/
    void updateClusterDependencies( unsigned int i_localClusterId );

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
     * @param i_memoryManager memory manager.
     * @param i_meshToClusters mapping from the mesh to the clusters.
     **/
    void addClusters( struct TimeStepping&               i_timeStepping,
                      struct MeshStructure*              i_meshStructure,
                      initializers::MemoryManager&       i_memoryManager );

    /**
     * Starts the communication thread.
     * Remark: This method has no effect when not compiled for communication thread support.
     **/
    void startCommunicationThread();

    /**
     * Stops the communication thread.
     * Remark: This method has no effect when not compiled for communication thread support.
     **/
    void stopCommunicationThread();

    /**
     * Advance in time until all clusters reach the next synchronization time.
     **/
    void advanceInTime( const double &i_synchronizationTime );

    /**
     * Gets the time tolerance of the time manager (1E-5 of the CFL time step width).
     **/
    double getTimeTolerance();
    
    /**
     * Distributes point sources pointers to clusters
     * 
     * @param i_pointSources Array of PointSources
     * @param i_numberOfLocalClusters Number of entries in i_pointSources
     */
    void setPointSourcesForClusters( sourceterm::ClusterMapping const* cms, sourceterm::PointSources const* pointSources );

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

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
    void pollForCommunication();

    static void* static_pollForCommunication(void* p) {
      static_cast<seissol::time_stepping::TimeManager*>(p)->pollForCommunication();
      return NULL;
    }
#endif

    void printComputationTime();
};

#endif

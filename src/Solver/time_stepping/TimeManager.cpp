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
 * Time Step width management in SeisSol.
 **/

#include "Parallel/MPI.h"

#include "TimeManager.h"
#include <Initializer/preProcessorMacros.fpp>
#include <Initializer/time_stepping/common.hpp>

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
#include <Parallel/Pin.h>
volatile bool g_executeCommThread;
volatile unsigned int* volatile g_handleRecvs;
volatile unsigned int* volatile g_handleSends;
pthread_t g_commThread;
#endif

seissol::time_stepping::TimeManager::TimeManager():
  m_logUpdates(std::numeric_limits<unsigned int>::max())
{
  m_loopStatistics.addRegion("computeLocalIntegration");
  m_loopStatistics.addRegion("computeNeighboringIntegration");
  m_loopStatistics.addRegion("computeDynamicRupture");
}

seissol::time_stepping::TimeManager::~TimeManager() {
  // free memory of the clusters
  for( unsigned l_cluster = 0; l_cluster < m_clusters.size(); l_cluster++ ) {
    delete m_clusters[l_cluster];
  }
}

void seissol::time_stepping::TimeManager::addClusters( struct TimeStepping&               i_timeStepping,
                                                       struct MeshStructure*              i_meshStructure,
                                                       initializers::MemoryManager&       i_memoryManager ) {
  SCOREP_USER_REGION( "addClusters", SCOREP_USER_REGION_TYPE_FUNCTION );

  // assert non-zero pointers
  assert( i_meshStructure         != NULL );

  // store the time stepping
  m_timeStepping = i_timeStepping;

  // iterate over local time clusters
  for( unsigned int l_cluster = 0; l_cluster < m_timeStepping.numberOfLocalClusters; l_cluster++ ) {
    struct MeshStructure          *l_meshStructure           = NULL;
    struct GlobalData             *l_globalData              = NULL;

    // get memory layout of this cluster
    i_memoryManager.getMemoryLayout( l_cluster,
                                     l_meshStructure,
                                     l_globalData
                                     );

    // add this time cluster
    m_clusters.push_back( new TimeCluster( l_cluster,
                                           m_timeStepping.clusterIds[l_cluster],
                                           l_meshStructure,
                                           l_globalData,
                                           &i_memoryManager.getLtsTree()->child(l_cluster),
                                           &i_memoryManager.getDynamicRuptureTree()->child(l_cluster),
                                           i_memoryManager.getLts(),
                                           i_memoryManager.getDynamicRupture(),
                                           &m_loopStatistics )
                        );
  }
}

void seissol::time_stepping::TimeManager::startCommunicationThread() {
#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
  g_executeCommThread = true;
  g_handleRecvs = (volatile unsigned int* volatile) malloc(sizeof(unsigned int) * m_timeStepping.numberOfLocalClusters);
  g_handleSends = (volatile unsigned int* volatile) malloc(sizeof(unsigned int) * m_timeStepping.numberOfLocalClusters);
  for ( unsigned int l_cluster = 0; l_cluster < m_timeStepping.numberOfLocalClusters; l_cluster++ ) {
    g_handleRecvs[l_cluster] = 0;
    g_handleSends[l_cluster] = 0;
  }
  pthread_create(&g_commThread, NULL, seissol::time_stepping::TimeManager::static_pollForCommunication, this);
#endif
}

void seissol::time_stepping::TimeManager::stopCommunicationThread() {
#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
  g_executeCommThread = false;
  pthread_join(g_commThread, NULL);
  free((void*)g_handleRecvs);
  free((void*)g_handleSends);
#endif
}

void seissol::time_stepping::TimeManager::updateClusterDependencies( unsigned int i_localClusterId ) {
  SCOREP_USER_REGION( "updateClusterDependencies", SCOREP_USER_REGION_TYPE_FUNCTION )

  // get time tolerance
  double l_timeTolerance = getTimeTolerance();

  // derive relevant clusters
  unsigned int l_lowerCluster = i_localClusterId;
  unsigned int l_upperCluster = i_localClusterId;

  if( i_localClusterId != 0 ) {
    l_lowerCluster--;
  }
  if( i_localClusterId < m_timeStepping.numberOfLocalClusters-1 ) {
    l_upperCluster++;
  }

  // iterate over the clusters
  for( unsigned int l_cluster = l_lowerCluster; l_cluster <= l_upperCluster; l_cluster++ ) {
    // get the relevant times
    double l_previousPredictionTime     = std::numeric_limits<double>::max();
    double l_previousFullUpdateTime     = std::numeric_limits<double>::max();
    double l_predictionTime             = m_clusters[l_cluster]->m_predictionTime;
    double l_fullUpdateTime             = m_clusters[l_cluster]->m_fullUpdateTime;
    double l_nextPredictionTime         = std::numeric_limits<double>::max();
    double l_nextUpcomingFullUpdateTime = std::numeric_limits<double>::max();

    if( l_cluster > 0 ) {
      l_previousPredictionTime = m_clusters[l_cluster-1]->m_predictionTime;
      l_previousFullUpdateTime = m_clusters[l_cluster-1]->m_fullUpdateTime;
    }

    if( l_cluster < m_timeStepping.numberOfLocalClusters - 1 ) {
      l_nextPredictionTime             = m_clusters[l_cluster+1]->m_predictionTime;
      l_nextUpcomingFullUpdateTime     = m_clusters[l_cluster+1]->m_fullUpdateTime + m_clusters[l_cluster+1]->timeStepWidth();
    }

    /*
     * Check if the cluster is eligible for a full update.
     *
     * Requirements for a full update:
     *  1) Cluster isn't queued already.
     *  2) Synchronization time isn't reached by now.
     *  3) Prediction time of the previous cluster is equal to the one of the current cluster.
     *  4) Prediction time of the cluster is equal to the desired time step width.
     *  5) Prediction time of the next cluster is greater or equal to the one of the current cluster.
     * _____________________________________________ _ _ _ _ _ _ _   _
     *    |             |             |             |             |   |
     * te | full update | full update | full update | prediction  |   |--- Status of the previous cluster
     * ___|_____________|_____________|_____________|_ _ _ _ _ _ _|  _|
     *
     * _________________ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _   _
     *                  |                                         |   |
     *  full update     | current prediction/planned full update  |   |--- Status of the current cluster.
     * _________________|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|  _|
     *
     * _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _   _
     *                                                            |   |
     *                          prediction                        |   |--- Status of the next cluster.
     * _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|  _|
     *
     */
     if( !m_clusters[l_cluster]->m_updatable.neighboringCopy && !m_clusters[l_cluster]->m_updatable.neighboringInterior                 &&  // 1)
          std::abs( l_fullUpdateTime - m_timeStepping.synchronizationTime )                        > l_timeTolerance                    &&  // 2)
          l_previousPredictionTime                                                                 > l_predictionTime - l_timeTolerance &&  // 3)
          std::abs( l_predictionTime - l_fullUpdateTime )                                          > l_timeTolerance                    &&  // 4)
          l_nextPredictionTime                                                                     > l_predictionTime - l_timeTolerance ) { // 5)
       // enqueue the cluster
#ifdef USE_MPI
       m_clusters[l_cluster]->m_updatable.neighboringCopy     = true;
       m_neighboringCopyQueue.push_back( m_clusters[l_cluster] );
#endif
       m_clusters[l_cluster]->m_updatable.neighboringInterior = true;
       m_neighboringInteriorQueue.push( m_clusters[l_cluster] );
     }

     /*
      * Check if the cluster is eligible for time prediction.
      *
      * Requirements for a prediction:
      *  1) Cluster isn't queued already.
      *  2) Synchronization time isn't reached by now.
      *  3) Previous cluster has reached the currents cluster current predicition time (doesn't require the current data anymore).
      *  4) Current cluster has used its own prediction data to perform a full update.
      *  5) Upcoming full update time of the next time cluster doesn't match the current prediction time <=> Next cluster has used the current clusters buffers to perform a full update.
      * _________________                                             _
      *    |             |                                             |
      * te | full update |                                             |--- Status of the previous cluster
      * ___|_____________|                                            _|
      *
      * _________________ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _   _
      *                  |                                         |   |
      *  full update     |             planned prediction          |   |--- Status of the current cluster.
      * _________________|_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|  _|
      *
      * _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _   _
      *                                                            |   |
      *                          prediction                        |   |--- Status of the next cluster.
      * _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|  _|
      */
      if( !m_clusters[l_cluster]->m_updatable.localCopy && !m_clusters[l_cluster]->m_updatable.localInterior                             &&  // 1)
          std::abs( l_fullUpdateTime - m_timeStepping.synchronizationTime )                        > l_timeTolerance                     &&  // 2)
          l_previousFullUpdateTime                                                                 > l_predictionTime - l_timeTolerance  &&  // 3)
          std::abs( l_fullUpdateTime - l_predictionTime )                                          < l_timeTolerance                     &&  // 4)
          l_nextUpcomingFullUpdateTime                                                             > l_predictionTime + l_timeTolerance  ) { // 5)
        // enqueue the cluster
#ifdef USE_MPI
        m_clusters[l_cluster]->m_updatable.localCopy = true;
        m_localCopyQueue.push_back( m_clusters[l_cluster] );
#endif
        m_clusters[l_cluster]->m_updatable.localInterior = true;
        m_localInteriorQueue.push( m_clusters[l_cluster] );

        // derive next time step width of the cluster
        unsigned int l_globalClusterId = m_timeStepping.clusterIds[l_cluster];
        // chop of at synchronization time
        m_clusters[l_cluster]->setTimeStepWidth( std::min( m_timeStepping.globalCflTimeStepWidths[l_globalClusterId],
                                                           m_timeStepping.synchronizationTime - m_clusters[l_cluster]->m_fullUpdateTime ) );

        // derive if the cluster is required to reset its lts buffers, reset sub time start and receive derivatives
        if( m_clusters[l_cluster]->m_numberOfFullUpdates % m_timeStepping.globalTimeStepRates[l_globalClusterId] == 0 ) {
          m_clusters[l_cluster]->m_resetLtsBuffers       = true;
          m_clusters[l_cluster]->m_subTimeStart          = 0;
        }
        else {
          m_clusters[l_cluster]->m_resetLtsBuffers       = false;
        }

#ifdef USE_MPI
        // TODO please check if this ifdef is correct

        // derive if the cluster is required to send its lts buffers
        if( (m_clusters[l_cluster]->m_numberOfFullUpdates + 1) % m_timeStepping.globalTimeStepRates[l_globalClusterId] == 0 ) {
          m_clusters[l_cluster]->m_sendLtsBuffers = true;
        }
        else {
          m_clusters[l_cluster]->m_sendLtsBuffers = false;
        }

        // derive if cluster is ready for synchronization
        if( std::abs( m_timeStepping.synchronizationTime - (m_clusters[l_cluster]->m_fullUpdateTime + m_clusters[l_cluster]->timeStepWidth()) ) < l_timeTolerance ) {
          m_clusters[l_cluster]->m_sendLtsBuffers = true;
        }
#endif // USE_MPI
      }
  }
}

void seissol::time_stepping::TimeManager::advanceInTime( const double &i_synchronizationTime ) {
  SCOREP_USER_REGION( "advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION )

  // asssert we are not moving back in time
  assert( m_timeStepping.synchronizationTime <= i_synchronizationTime );

  m_timeStepping.synchronizationTime = i_synchronizationTime;

  // iterate over all clusters and set default values
  for( unsigned int l_cluster = 0; l_cluster < m_timeStepping.numberOfLocalClusters; l_cluster++ ) {
#ifdef USE_MPI
    m_clusters[l_cluster]->m_updatable.localCopy           = false;
    m_clusters[l_cluster]->m_updatable.neighboringCopy     = false;
#endif
    m_clusters[l_cluster]->m_updatable.localInterior       = false;
    m_clusters[l_cluster]->m_updatable.neighboringInterior = false;

    m_clusters[l_cluster]->m_resetLtsBuffers               = true;
    m_clusters[l_cluster]->setTimeStepWidth(0.);
    m_clusters[l_cluster]->m_subTimeStart                  = 0;
    m_clusters[l_cluster]->m_numberOfFullUpdates           = 0;
  }

  // initialize prediction queues
  for( unsigned int l_cluster = 0; l_cluster < m_timeStepping.numberOfLocalClusters; l_cluster++ ) {
    updateClusterDependencies(l_cluster);
  }

  // iterate until all queues are empty and the next synchronization point in time is reached
  while( !( m_localCopyQueue.empty()       && m_localInteriorQueue.empty() &&
            m_neighboringCopyQueue.empty() && m_neighboringInteriorQueue.empty() ) ) {
#ifdef USE_MPI
    // iterate over all items of the local copy queue and update everything possible
    for( std::list<TimeCluster*>::iterator l_cluster = m_localCopyQueue.begin(); l_cluster != m_localCopyQueue.end(); ) {
      if( (*l_cluster)->computeLocalCopy() ) {
        unsigned int l_clusterId = (*l_cluster)->m_clusterId;
        l_cluster = m_localCopyQueue.erase( l_cluster );
        updateClusterDependencies( l_clusterId );
      }
      else l_cluster++;
    }

    // iterate over all items of the neighboring copy queue and update everything possible
    for( std::list<TimeCluster*>::iterator l_cluster = m_neighboringCopyQueue.begin(); l_cluster != m_neighboringCopyQueue.end(); ) {
      if( (*l_cluster)->computeNeighboringCopy() ) {
        unsigned int l_clusterId = (*l_cluster)->m_clusterId;
        l_cluster = m_neighboringCopyQueue.erase( l_cluster );
        updateClusterDependencies( l_clusterId );
      }
      else l_cluster++;
    }
#endif

    // update a single interior region (if present) with local updates
    if( !m_localInteriorQueue.empty() ) {
      TimeCluster *l_timeCluster = m_localInteriorQueue.top();
      l_timeCluster->computeLocalInterior();
      m_localInteriorQueue.pop();
      updateClusterDependencies(l_timeCluster->m_clusterId);
    }

    // update a single interior region (if present) with neighboring updates
    if( !m_neighboringInteriorQueue.empty() ) {
      TimeCluster *l_timeCluster = m_neighboringInteriorQueue.top();
      l_timeCluster->computeNeighboringInterior();
      m_neighboringInteriorQueue.pop();
      updateClusterDependencies(l_timeCluster->m_clusterId);
    }

    // print progress of largest time cluster
    if( m_clusters[m_timeStepping.numberOfLocalClusters-1]->m_numberOfFullUpdates != m_logUpdates &&
        m_clusters[m_timeStepping.numberOfLocalClusters-1]->m_numberOfFullUpdates % 100 == 0 ) {
      m_logUpdates = m_clusters[m_timeStepping.numberOfLocalClusters-1]->m_numberOfFullUpdates;

      const int rank = MPI::mpi.rank();

      logInfo(rank) << "#max-updates since sync: " << m_logUpdates
                         << " @ "                  << m_clusters[m_timeStepping.numberOfLocalClusters-1]->m_fullUpdateTime;
    }
  }
}

void seissol::time_stepping::TimeManager::printComputationTime()
{
#ifdef USE_MPI
  m_loopStatistics.printSummary(MPI::mpi.comm());
#endif
  m_loopStatistics.writeSamples();
}

double seissol::time_stepping::TimeManager::getTimeTolerance() {
  return 1E-5 * m_timeStepping.globalCflTimeStepWidths[0];
}

void seissol::time_stepping::TimeManager::setPointSourcesForClusters( sourceterm::ClusterMapping const* cms, sourceterm::PointSources const* pointSources )
{
  for (unsigned cluster = 0; cluster < m_clusters.size(); ++cluster) {
    m_clusters[cluster]->setPointSources( cms[cluster].cellToSources,
                                          cms[cluster].numberOfMappings,
                                          &pointSources[cluster] );
  }
}

void seissol::time_stepping::TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter)
{
  for (unsigned cluster = 0; cluster < m_clusters.size(); ++cluster) {
    m_clusters[cluster]->setReceiverCluster(receiverWriter.receiverCluster(cluster));
  }
}

void seissol::time_stepping::TimeManager::setInitialTimes( double i_time ) {
  assert( i_time >= 0 );

  for( unsigned int l_cluster = 0; l_cluster < m_clusters.size(); l_cluster++ ) {
    m_clusters[l_cluster]->m_predictionTime = i_time;
    m_clusters[l_cluster]->m_fullUpdateTime = i_time;
    m_clusters[l_cluster]->m_receiverTime   = i_time;
  }
}

void seissol::time_stepping::TimeManager::setTv(double tv) {
  for( unsigned int l_cluster = 0; l_cluster < m_clusters.size(); l_cluster++ ) {
    m_clusters[l_cluster]->setTv(tv);
  }
}

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
void seissol::time_stepping::TimeManager::pollForCommunication() {
  // pin this thread to the last core
  volatile unsigned int l_signalSum = 0;

  parallel::pinToFreeCPUs();

  //logInfo(0) << "Launching communication thread on OS core id:" << l_numberOfHWThreads;

  // now let's enter the polling loop
  while (g_executeCommThread == true || l_signalSum > 0) {
    for( unsigned int l_cluster = 0; l_cluster < m_clusters.size(); l_cluster++ ) {
      if (g_handleRecvs[l_cluster] == 1) {
        m_clusters[l_cluster]->startReceiveGhostLayer();
        g_handleRecvs[l_cluster] = 2;
      }
      if (g_handleSends[l_cluster] == 1) {
        m_clusters[l_cluster]->startSendCopyLayer();
        g_handleSends[l_cluster] = 2;
      }
      if (g_handleRecvs[l_cluster] == 2) {
        m_clusters[l_cluster]->pollForGhostLayerReceives();
      }
      if (g_handleSends[l_cluster] == 2) {
        m_clusters[l_cluster]->pollForCopyLayerSends();
      }
    }
    l_signalSum = 0;
    for( unsigned int l_cluster = 0; l_cluster < m_clusters.size(); l_cluster++ ) {
      l_signalSum += g_handleRecvs[l_cluster];
      l_signalSum += g_handleSends[l_cluster];
    }
  }
}
#endif


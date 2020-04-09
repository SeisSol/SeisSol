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

#include <numeric>
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
  for(unsigned l_cluster = 0; l_cluster < clusters.size(); l_cluster++ ) {
    delete clusters[l_cluster];
  }
}

void seissol::time_stepping::TimeManager::addClusters(struct TimeStepping& i_timeStepping,
                                                      struct MeshStructure* i_meshStructure,
                                                      initializers::MemoryManager& i_memoryManager) {
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

    const unsigned int l_globalClusterId = m_timeStepping.clusterIds[l_cluster];
    // chop of at synchronization time
    const auto timeStepSize = m_timeStepping.globalCflTimeStepWidths[l_globalClusterId];
    const auto timeStepRate =  static_cast<int>(m_timeStepping.globalTimeStepRates[l_globalClusterId]);
        // add this time cluster
    const auto layerTypes = {Copy, Interior};
    for (auto type : layerTypes) {
      clusters.push_back(new TimeCluster(
          l_cluster,
          m_timeStepping.clusterIds[l_cluster],
          timeStepSize,
          timeStepRate,
          getTimeTolerance(),
          l_globalData,
          &i_memoryManager.getLtsTree()->child(l_cluster).child(type),
          &i_memoryManager.getDynamicRuptureTree()->child(l_cluster).child(type),
          i_memoryManager.getLts(),
          i_memoryManager.getDynamicRupture(),
          &m_loopStatistics)
      );
    }
    auto* interior = clusters[clusters.size() - 1];
    auto* copy = clusters[clusters.size() - 2];
    // Copy/interior with same timestep are neighbors
    interior->connect(*copy);

    // Connect new copy/interior to previous two copy/interior
    // Then all clusters that are neighboring are connected.
    // Note: Only clusters with a distance of 1 time step factor
    // are connected.
    if (l_cluster > 0) {
      assert(clusters.size() >= 4);
      for (int i = 0; i < 2; ++i)  {
        copy->connect(
            *clusters[clusters.size() - 2 - i - 1]
        );
        interior->connect(
            *clusters[clusters.size() - 2 - i - 1]
        );
      }
    }

#ifdef USE_MPI
    // Create ghost time clusters for MPI
    const int globalClusterId = static_cast<int>(m_timeStepping.clusterIds[l_cluster]);
    /*for (int otherGlobalClusterId = std::max(globalClusterId - 1, 0);
         otherGlobalClusterId < std::min(globalClusterId + 2, static_cast<int>(m_timeStepping.numberOfGlobalClusters));
         ++otherGlobalClusterId) {
         */
    for (unsigned int otherGlobalClusterId = 0; otherGlobalClusterId < m_timeStepping.numberOfGlobalClusters; ++otherGlobalClusterId) {
      const bool hasNeighborRegions = std::any_of(l_meshStructure->neighboringClusters,
      l_meshStructure->neighboringClusters + l_meshStructure->numberOfRegions,
      [otherGlobalClusterId](const auto& neighbor) {
        return neighbor[1] == otherGlobalClusterId;
      });
      if (hasNeighborRegions) {
        const auto otherTimeStepSize = m_timeStepping.globalCflTimeStepWidths[otherGlobalClusterId];
        ghostClusters.push_back(
          std::make_unique<GhostTimeCluster>(
              otherTimeStepSize,
              timeStepRate,
              getTimeTolerance(),
              otherGlobalClusterId,
              l_meshStructure)
        );
        // Connect with previous copy layer.
        ghostClusters.back()->connect(*clusters[clusters.size() - 2]);
      }
    }
#endif
  }
}

void seissol::time_stepping::TimeManager::startCommunicationThread() {
  /*
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
   */
}

void seissol::time_stepping::TimeManager::stopCommunicationThread() {
  /*
#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
  g_executeCommThread = false;
  pthread_join(g_commThread, NULL);
  free((void*)g_handleRecvs);
  free((void*)g_handleSends);
#endif
   */
}

void seissol::time_stepping::TimeManager::updateClusterDependencies( unsigned int i_localClusterId ) {
  /*
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
    double l_predictionTime             = clusters[l_cluster]->m_predictionTime;
    double l_fullUpdateTime             = clusters[l_cluster]->m_fullUpdateTime;
    double l_nextPredictionTime         = std::numeric_limits<double>::max();
    double l_nextUpcomingFullUpdateTime = std::numeric_limits<double>::max();

    if( l_cluster > 0 ) {
      l_previousPredictionTime = clusters[l_cluster-1]->m_predictionTime;
      l_previousFullUpdateTime = clusters[l_cluster-1]->m_fullUpdateTime;
    }

    if( l_cluster < m_timeStepping.numberOfLocalClusters - 1 ) {
      l_nextPredictionTime             = clusters[l_cluster+1]->m_predictionTime;
      l_nextUpcomingFullUpdateTime     = clusters[l_cluster+1]->m_fullUpdateTime + clusters[l_cluster+1]->timeStepWidth();
    }*/

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
     /*if( !clusters[l_cluster]->m_updatable.neighboringCopy && !clusters[l_cluster]->m_updatable.neighboringInterior                 &&  // 1)
          std::abs( l_fullUpdateTime - m_timeStepping.synchronizationTime )                        > l_timeTolerance                    &&  // 2)
          l_previousPredictionTime                                                                 > l_predictionTime - l_timeTolerance &&  // 3)
          std::abs( l_predictionTime - l_fullUpdateTime )                                          > l_timeTolerance                    &&  // 4)
          l_nextPredictionTime                                                                     > l_predictionTime - l_timeTolerance ) { // 5)
       // enqueue the cluster
#ifdef USE_MPI
       clusters[l_cluster]->m_updatable.neighboringCopy     = true;
       m_neighboringCopyQueue.push_back( clusters[l_cluster] );
#endif
       clusters[l_cluster]->m_updatable.neighboringInterior = true;
       m_neighboringInteriorQueue.push( clusters[l_cluster] );
     }*/

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
      /*if( !clusters[l_cluster]->m_updatable.localCopy && !clusters[l_cluster]->m_updatable.localInterior                             &&  // 1)
          std::abs( l_fullUpdateTime - m_timeStepping.synchronizationTime )                        > l_timeTolerance                     &&  // 2)
          l_previousFullUpdateTime                                                                 > l_predictionTime - l_timeTolerance  &&  // 3)
          std::abs( l_fullUpdateTime - l_predictionTime )                                          < l_timeTolerance                     &&  // 4)
          l_nextUpcomingFullUpdateTime                                                             > l_predictionTime + l_timeTolerance  ) { // 5)
        // enqueue the cluster
#ifdef USE_MPI
        clusters[l_cluster]->m_updatable.localCopy = true;
        m_localCopyQueue.push_back( clusters[l_cluster] );
#endif
        clusters[l_cluster]->m_updatable.localInterior = true;
        m_localInteriorQueue.push( clusters[l_cluster] );

        // derive next time step width of the cluster
        unsigned int l_globalClusterId = m_timeStepping.clusterIds[l_cluster];
        // chop of at synchronization time
        clusters[l_cluster]->setTimeStepWidth( std::min( m_timeStepping.globalCflTimeStepWidths[l_globalClusterId],
                                                           m_timeStepping.synchronizationTime - clusters[l_cluster]->m_fullUpdateTime ) );

        // derive if the cluster is required to reset its lts buffers, reset sub time start and receive derivatives
        if( clusters[l_cluster]->m_numberOfFullUpdates % m_timeStepping.globalTimeStepRates[l_globalClusterId] == 0 ) {
          clusters[l_cluster]->m_resetLtsBuffers       = true;
          clusters[l_cluster]->m_subTimeStart          = 0;
        }
        else {
          clusters[l_cluster]->m_resetLtsBuffers       = false;
        }

#ifdef USE_MPI
        // TODO please check if this ifdef is correct

        // derive if the cluster is required to send its lts buffers
        if( (clusters[l_cluster]->m_numberOfFullUpdates + 1) % m_timeStepping.globalTimeStepRates[l_globalClusterId] == 0 ) {
          clusters[l_cluster]->m_sendLtsBuffers = true;
        }
        else {
          clusters[l_cluster]->m_sendLtsBuffers = false;
        }

        // derive if cluster is ready for synchronization
        if( std::abs( m_timeStepping.synchronizationTime - (clusters[l_cluster]->m_fullUpdateTime + clusters[l_cluster]->timeStepWidth()) ) < l_timeTolerance ) {
          clusters[l_cluster]->m_sendLtsBuffers = true;
        }
#endif // USE_MPI
      }
  }*/
}

void seissol::time_stepping::TimeManager::advanceInTime(const double &synchronizationTime) {
  SCOREP_USER_REGION( "advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION )

  // We should always move forward in time
  assert(m_timeStepping.synchronizationTime <= synchronizationTime);

  m_timeStepping.synchronizationTime = synchronizationTime;
  std::cout << seissol::MPI::mpi.rank() << " new sync time = " << synchronizationTime << std::endl;

  // TODO(Lukas) WTF...
  //MPI_Comm newComm;
  //MPI_Comm_dup(seissol::MPI::mpi.comm(), &newComm);
  //seissol::MPI::mpi.setComm(newComm);

  for (auto* cluster : clusters) {
    cluster->reset();
    cluster->updateSyncTime(synchronizationTime);
  }
  for (auto& ghostCluster : ghostClusters) {
    // TODO(Lukas) Not sure about cancel + reset
    ghostCluster->cancelPendingMessages();
    ghostCluster->reset();

    ghostCluster->updateSyncTime(synchronizationTime);
  }
  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());

  bool finished = false; // Is true, once all clusters reached next sync point
  while (!finished) {
    finished = true;
    for (auto* cluster : clusters) {
      // TODO(Lukas) Remove
      cluster->resetBuffersOld = cluster->numberOfTimeSteps % m_timeStepping.globalTimeStepRates[cluster->m_globalClusterId] == 0;

      // A cluster yields once it is blocked by other cluster.
      bool yield = false;
      do {
        yield = cluster->act();
        // Check ghost cells often for communication progress
        // Note: This replaces the need for a communication thread.
        for (auto& ghostCluster : ghostClusters) {
          ghostCluster->act();
          finished = finished && ghostCluster->synced();
        }
      if (!(yield || cluster->synced())) {
              //std::cout << seissol::MPI::mpi.rank() <<  "Ghost cluster act" << std::endl;
          }
      } while (!(yield || cluster->synced()));
      finished = finished && cluster->synced();
    }
  }

  // iterate over all clusters and set default values
  /*for( unsigned int l_cluster = 0; l_cluster < m_timeStepping.numberOfLocalClusters; l_cluster++ ) {
#ifdef USE_MPI
    clusters[l_cluster]->m_updatable.localCopy           = false;
    clusters[l_cluster]->m_updatable.neighboringCopy     = false;
#endif
    clusters[l_cluster]->m_updatable.localInterior       = false;
    clusters[l_cluster]->m_updatable.neighboringInterior = false;

    clusters[l_cluster]->m_resetLtsBuffers               = true;
    clusters[l_cluster]->setTimeStepWidth(0.);
    clusters[l_cluster]->m_subTimeStart                  = 0;
    clusters[l_cluster]->m_numberOfFullUpdates           = 0;
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
    if( clusters[m_timeStepping.numberOfLocalClusters-1]->m_numberOfFullUpdates != m_logUpdates &&
        clusters[m_timeStepping.numberOfLocalClusters-1]->m_numberOfFullUpdates % 100 == 0 ) {
      m_logUpdates = clusters[m_timeStepping.numberOfLocalClusters-1]->m_numberOfFullUpdates;

      const int rank = MPI::mpi.rank();

      logInfo(rank) << "#max-updates since sync: " << m_logUpdates
                         << " @ "                  << clusters[m_timeStepping.numberOfLocalClusters-1]->m_fullUpdateTime;
    }
  }*/
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
  for (unsigned cluster = 0; cluster < clusters.size(); ++cluster) {
    clusters[cluster]->setPointSources(cms[cluster].cellToSources,
                                       cms[cluster].numberOfMappings,
                                       &pointSources[cluster] );
  }
}

void seissol::time_stepping::TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter)
{
  for (unsigned cluster = 0; cluster < clusters.size(); ++cluster) {
    clusters[cluster]->setReceiverCluster(receiverWriter.receiverCluster(cluster));
  }
}

void seissol::time_stepping::TimeManager::setInitialTimes( double i_time ) {
  assert( i_time >= 0 );

  for(unsigned int l_cluster = 0; l_cluster < clusters.size(); l_cluster++ ) {
    // TODO set initial times for checkpointing
    //clusters[l_cluster]->m_predictionTime = i_time;
    //clusters[l_cluster]->m_fullUpdateTime = i_time;
    clusters[l_cluster]->m_receiverTime   = i_time;
  }
}

void seissol::time_stepping::TimeManager::setTv(double tv) {
  for(auto& cluster : clusters) {
    cluster->setTv(tv);
  }
}
void seissol::time_stepping::TimeManager::cancelPendingMessages() {
  for (auto& cluster : ghostClusters) {
    cluster->cancelPendingMessages();
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
    for( unsigned int l_cluster = 0; l_cluster < clusters.size(); l_cluster++ ) {
      if (g_handleRecvs[l_cluster] == 1) {
        clusters[l_cluster]->startReceiveGhostLayer();
        g_handleRecvs[l_cluster] = 2;
      }
      if (g_handleSends[l_cluster] == 1) {
        clusters[l_cluster]->startSendCopyLayer();
        g_handleSends[l_cluster] = 2;
      }
      if (g_handleRecvs[l_cluster] == 2) {
        clusters[l_cluster]->pollForGhostLayerReceives();
      }
      if (g_handleSends[l_cluster] == 2) {
        clusters[l_cluster]->pollForCopyLayerSends();
      }
    }
    l_signalSum = 0;
    for( unsigned int l_cluster = 0; l_cluster < clusters.size(); l_cluster++ ) {
      l_signalSum += g_handleRecvs[l_cluster];
      l_signalSum += g_handleSends[l_cluster];
    }
  }
}
#endif


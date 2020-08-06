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
  for (unsigned int localClusterId = 0; localClusterId < m_timeStepping.numberOfLocalClusters; localClusterId++) {
    struct MeshStructure *l_meshStructure = nullptr;
    struct GlobalData *l_globalData = nullptr;
    i_memoryManager.getMemoryLayout(localClusterId,
                                    l_meshStructure,
                                    l_globalData
                                     );

    const unsigned int l_globalClusterId = m_timeStepping.clusterIds[localClusterId];
    // chop of at synchronization time
    const auto timeStepSize = m_timeStepping.globalCflTimeStepWidths[l_globalClusterId];
    const long timeStepRate = ipow(static_cast<long>(m_timeStepping.globalTimeStepRates[0]),
         static_cast<long>(l_globalClusterId));
    const auto layerTypes = {Copy, Interior};
    for (auto type : layerTypes) {
      // We print progress only if it is the cluster with the largest time step on each rank.
      // This does not mean that it is the largest cluster globally!
      const bool printProgress = (localClusterId == 0) && (type == Interior);
      clusters.push_back(new TimeCluster(
          localClusterId,
          m_timeStepping.clusterIds[localClusterId],
          timeStepSize,
          timeStepRate,
          getTimeTolerance(),
          printProgress,
          l_globalData,
          &i_memoryManager.getLtsTree()->child(localClusterId).child(type),
          &i_memoryManager.getDynamicRuptureTree()->child(localClusterId).child(type),
          i_memoryManager.getLts(),
          i_memoryManager.getDynamicRupture(),
          &m_loopStatistics)
      );
    }
    auto* interior = clusters[clusters.size() - 1];
    auto* copy = clusters[clusters.size() - 2];

    // Mark copy layers as higher priority layers.
    constexpr int priorityLow = -1;
    constexpr int priorityHigh = 42;
    interior->setPriority(priorityLow);
    copy->setPriority(priorityHigh);

    // Copy/interior with same timestep are neighbors
    interior->connect(*copy);

    // Connect new copy/interior to previous two copy/interior
    // Then all clusters that are neighboring are connected.
    // Note: Only clusters with a distance of 1 time step factor
    // are connected.
    if (localClusterId > 0) {
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
    const int globalClusterId = static_cast<int>(m_timeStepping.clusterIds[localClusterId]);
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
          assert(otherGlobalClusterId >= std::max(globalClusterId - 1, 0));
          assert(otherGlobalClusterId < std::min(globalClusterId +2, static_cast<int>(m_timeStepping.numberOfGlobalClusters)));
        const auto otherTimeStepSize = m_timeStepping.globalCflTimeStepWidths[otherGlobalClusterId];
        const auto otherTimeStepRate = ipow(2l, static_cast<long>(otherGlobalClusterId));

        // TODO(Lukas) Should also pass own timeStepRate for checking whether to send etc
        ghostClusters.push_back(
          std::make_unique<GhostTimeCluster>(
              otherTimeStepSize,
              otherTimeStepRate,
              getTimeTolerance(),
              globalClusterId,
              otherGlobalClusterId,
              l_meshStructure)
        );
        // Connect with previous copy layer.
        ghostClusters.back()->connect(*clusters[clusters.size() - 2]);
      }
    }
#endif
  }

  // Sort clusters by priority (higher ones first).
  auto prioritySorter = [](const auto& a, const auto& b) {
    return a->getPriority() < b->getPriority();
  };
  std::sort(clusters.begin(), clusters.end(), prioritySorter);
  std::sort(ghostClusters.begin(), ghostClusters.end(), prioritySorter);
}



void seissol::time_stepping::TimeManager::advanceInTime(const double &synchronizationTime) {
  SCOREP_USER_REGION( "advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION )

  // We should always move forward in time
  assert(m_timeStepping.synchronizationTime <= synchronizationTime);

  m_timeStepping.synchronizationTime = synchronizationTime;
  std::cout << seissol::MPI::mpi.rank() << " new sync time = " << synchronizationTime << std::endl;

  for (auto* cluster : clusters) {
    cluster->updateSyncTime(synchronizationTime);
    cluster->reset();
  }
  // TODO(Lukas) Remove.
  //seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());

  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->updateSyncTime(synchronizationTime);
    ghostCluster->reset();
  }
  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());

  bool finished = false; // Is true, once all clusters reached next sync point
  while (!finished) {
    finished = true;
    for (auto* cluster : clusters) {
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
      } while (!(yield || cluster->synced()));
      finished = finished && cluster->synced();
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


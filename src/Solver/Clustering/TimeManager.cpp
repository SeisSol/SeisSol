/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#include "Parallel/Helper.hpp"
#include "ResultWriter/ClusteringWriter.h"
#include "SeisSol.h"
#include "Solver/Clustering/Communication/CommunicationManager.h"
#include "TimeManager.h"
#include <AbstractAPI.h>
#include <Common/Executor.hpp>
#include <DynamicRupture/Output/OutputManager.hpp>
#include <Initializer/MemoryManager.h>
#include <Initializer/tree/Layer.hpp>
#include <Initializer/typedefs.hpp>
#include <Kernels/PointSourceCluster.h>
#include <ResultWriter/ReceiverWriter.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/AbstractGhostTimeCluster.h>
#include <Solver/Clustering/Communication/GhostTimeClusterFactory.h>
#include <Solver/Clustering/Computation/DynamicRuptureCluster.hpp>
#include <Solver/Clustering/Computation/TimeCluster.h>
#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>
#include <xdmfwriter/scorep_wrapper.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::time_stepping {

TimeManager::TimeManager(seissol::SeisSol& seissolInstance)
    : logUpdates(std::numeric_limits<unsigned int>::max()), seissolInstance(seissolInstance),
      actorStateStatisticsManager(loopStatistics) {
  loopStatistics.addRegion("computeLocalIntegration");
  loopStatistics.addRegion("computeNeighboringIntegration");
  loopStatistics.addRegion("computeDynamicRupture");
  loopStatistics.addRegion("computePointSources");

  loopStatistics.enableSampleOutput(
      seissolInstance.getSeisSolParameters().output.loopStatisticsNetcdfOutput);
}

void TimeManager::addClusters(TimeStepping& timeStepping,
                              MeshStructure* meshStructure,
                              initializer::MemoryManager& memoryManager,
                              bool usePlasticity) {
  SCOREP_USER_REGION("addClusters", SCOREP_USER_REGION_TYPE_FUNCTION);
  std::vector<std::unique_ptr<AbstractTimeCluster>> ghostClusters;
  // assert non-zero pointers
  assert(meshStructure != NULL);

  // store the time stepping
  this->timeStepping = timeStepping;

  auto clusteringWriter = writer::ClusteringWriter(memoryManager.getOutputPrefix());

  bool foundDynamicRuptureCluster = false;

  // iterate over local time clusters
  for (unsigned int localClusterId = 0; localClusterId < timeStepping.numberOfLocalClusters;
       localClusterId++) {
    // get memory layout of this cluster
    auto [meshStructure, globalData] = memoryManager.getMemoryLayout(localClusterId);

    const unsigned int globalClusterId = timeStepping.clusterIds[localClusterId];
    // chop off at synchronization time
    const auto timeStepSize = timeStepping.globalCflTimeStepWidths[globalClusterId];
    const long timeStepRate = ipow(static_cast<long>(timeStepping.globalTimeStepRates[0]),
                                   static_cast<long>(globalClusterId));

    // Dynamic rupture
    auto& dynRupTree = memoryManager.getDynamicRuptureTree()->child(localClusterId);
    // Note: We need to include the Ghost part, as we need to compute its DR part as well.
    const long numberOfDynRupCells = dynRupTree.child(Interior).getNumberOfCells() +
                                     dynRupTree.child(Copy).getNumberOfCells() +
                                     dynRupTree.child(Ghost).getNumberOfCells();

    bool isFirstDynamicRuptureCluster = false;
    if (!foundDynamicRuptureCluster && numberOfDynRupCells > 0) {
      foundDynamicRuptureCluster = true;
      isFirstDynamicRuptureCluster = true;
    }

    for (auto type : {Copy, Interior}) {
      const auto offsetMonitoring = type == Interior ? 0 : timeStepping.numberOfGlobalClusters;
      // We print progress only if it is the cluster with the largest time step on each rank.
      // This does not mean that it is the largest cluster globally!
      const bool printProgress =
          (localClusterId == timeStepping.numberOfLocalClusters - 1) && (type == Interior);
      const auto profilingId = globalClusterId + offsetMonitoring;
      auto* layerData = &memoryManager.getLtsTree()->child(localClusterId).child(type);
      auto* dynRupData = &dynRupTree.child(type);
      clusters.push_back(
          std::make_unique<TimeCluster>(localClusterId,
                                        globalClusterId,
                                        profilingId,
                                        usePlasticity,
                                        type,
                                        timeStepSize,
                                        timeStepRate,
                                        printProgress,
                                        globalData,
                                        layerData,
                                        memoryManager.getLts(),
                                        seissolInstance,
                                        &loopStatistics,
                                        &actorStateStatisticsManager.addCluster(profilingId)));
      clustersDR.push_back(std::make_unique<DynamicRuptureCluster>(
          timeStepSize,
          timeStepRate,
          profilingId,
          dynRupData,
          memoryManager.getDynamicRupture(),
          globalData,
          memoryManager.getFrictionLaw(),
          memoryManager.getFrictionLawDevice(),
          memoryManager.getFaultOutputManager(),
          seissolInstance,
          &loopStatistics,
          &actorStateStatisticsManager.addCluster(profilingId)));

      const auto clusterSize = layerData->getNumberOfCells();
      const auto dynRupSize = dynRupData->getNumberOfCells();
      // Add writer to output
      clusteringWriter.addCluster(profilingId, localClusterId, type, clusterSize, dynRupSize);
    }
    auto& interior = clusters[clusters.size() - 1];
    auto& copy = clusters[clusters.size() - 2];
    auto& interiorDR = clustersDR[clustersDR.size() - 1];
    auto& copyDR = clustersDR[clustersDR.size() - 2];

    // Mark copy layers as higher priority layers.
    interior->setPriority(ActorPriority::Low);
    copy->setPriority(ActorPriority::High);

    interiorDR->setPriority(ActorPriority::Low);
    copyDR->setPriority(ActorPriority::Low);

    // Copy/interior with same timestep are neighbors
    interior->connect(*copy);

    // connect clusters to DR clusters
    interior->connect(*interiorDR);
    copy->connect(*interiorDR);
    copy->connect(*copyDR);

    // Connect new copy/interior to previous two copy/interior
    // Then all clusters that are neighboring are connected.
    // Note: Only clusters with a distance of 1 time step factor
    // are connected.
    if (localClusterId > 0) {
      assert(clusters.size() >= 4);
      for (int i = 0; i < 2; ++i) {
        copy->connect(*clusters[clusters.size() - 2 - i - 1]);
        interior->connect(*clusters[clusters.size() - 2 - i - 1]);
      }
    }

    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    auto& [increaseManager, decreaseManager] = timeMirrorManagers;

    increaseManager.setTimeClusterVector(&clusters);
    decreaseManager.setTimeClusterVector(&clusters);
#ifdef USE_MPI
    // Create ghost time clusters for MPI
    const auto preferredDataTransferMode = MPI::mpi.getPreferredDataTransferMode();
    const auto persistent = usePersistentMpi();
    GhostTimeClusterFactory::prepare(preferredDataTransferMode,
                                     timeStepping.numberOfGlobalClusters);

    for (unsigned int otherGlobalClusterId = 0;
         otherGlobalClusterId < timeStepping.numberOfGlobalClusters;
         ++otherGlobalClusterId) {
      const bool hasNeighborRegions =
          std::any_of(meshStructure->neighboringClusters,
                      meshStructure->neighboringClusters + meshStructure->numberOfRegions,
                      [otherGlobalClusterId](const auto& neighbor) {
                        return static_cast<unsigned>(neighbor[1]) == otherGlobalClusterId;
                      });
      if (hasNeighborRegions) {
        assert(otherGlobalClusterId >= std::max(static_cast<int>(globalClusterId - 1), 0));
        assert(otherGlobalClusterId <
               std::min(globalClusterId + 2, timeStepping.numberOfGlobalClusters));
        const auto otherTimeStepSize = timeStepping.globalCflTimeStepWidths[otherGlobalClusterId];
        const long otherTimeStepRate = ipow(static_cast<long>(timeStepping.globalTimeStepRates[0]),
                                            static_cast<long>(otherGlobalClusterId));

        auto ghostCluster = GhostTimeClusterFactory::get(otherTimeStepSize,
                                                         otherTimeStepRate,
                                                         globalClusterId,
                                                         otherGlobalClusterId,
                                                         meshStructure,
                                                         preferredDataTransferMode,
                                                         persistent);
        ghostClusters.push_back(std::move(ghostCluster));

        // Connect with previous copy layer.
        ghostClusters.back()->connect(*copy);
      }
    }
#endif
  }

  clusteringWriter.write();

  // Sort clusters by time step size in increasing order
  auto rateSorter = [](const auto& a, const auto& b) {
    return a->getTimeStepRate() < b->getTimeStepRate();
  };
  std::sort(clusters.begin(), clusters.end(), rateSorter);
  std::sort(clustersDR.begin(), clustersDR.end(), rateSorter);

  for (const auto& cluster : clusters) {
    if (cluster->getPriority() == ActorPriority::High) {
      highPrioClusters.emplace_back(cluster.get());
    } else {
      lowPrioClusters.emplace_back(cluster.get());
    }
  }
  for (const auto& cluster : clustersDR) {
    if (cluster->getPriority() == ActorPriority::High) {
      highPrioClustersDR.emplace_back(cluster.get());
    } else {
      lowPrioClustersDR.emplace_back(cluster.get());
    }
  }

  std::sort(ghostClusters.begin(), ghostClusters.end(), rateSorter);

  if (seissol::useCommThread(MPI::mpi)) {
    communicationManager = std::make_unique<ThreadedCommunicationManager>(
        std::move(ghostClusters), &seissolInstance.getPinning());
  } else {
    communicationManager = std::make_unique<SerialCommunicationManager>(std::move(ghostClusters));
  }

  auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
  auto& [increaseManager, decreaseManager] = timeMirrorManagers;

  auto ghostClusterPointer = communicationManager->getGhostClusters();

  increaseManager.setGhostClusterVector(ghostClusterPointer);
  decreaseManager.setGhostClusterVector(ghostClusterPointer);
}

void TimeManager::setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
  this->faultOutputManager = faultOutputManager;
  for (auto& cluster : clustersDR) {
    cluster->setFaultOutputManager(faultOutputManager);
  }
}

seissol::dr::output::OutputManager* TimeManager::getFaultOutputManager() {
  assert(faultOutputManager != nullptr);
  return faultOutputManager;
}

void TimeManager::advanceInTime(const double& synchronizationTime) {
  SCOREP_USER_REGION("advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION)

  // We should always move forward in time
  assert(timeStepping.synchronizationTime <= synchronizationTime);

  timeStepping.synchronizationTime = synchronizationTime;

  for (auto& cluster : clusters) {
    cluster->setSyncTime(synchronizationTime);
    cluster->reset();
  }
  for (auto& cluster : clustersDR) {
    cluster->setSyncTime(synchronizationTime);
    cluster->reset();
  }

  communicationManager->reset(synchronizationTime);

  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->putProfilingMark("advanceInTime", device::ProfilingColors::Blue);
#endif

  const auto phaseExecute = [&]() {
    // Move all clusters from RestartAfterSync to Corrected
    // Does not involve any computations
    for (auto& cluster : clusters) {
      cluster->act();
    }

    bool finished = false; // Is true, once all clusters reached next sync point
    while (!finished) {
      finished = true;
      communicationManager->progression();

      // Update all high priority clusters
      for (int i = 0; i < highPrioClusters.size(); ++i) {
        std::for_each(highPrioClusters.begin(), highPrioClusters.end(), [&](auto& cluster) {
          communicationManager->progression();
          cluster->act();
        });
        std::for_each(highPrioClustersDR.begin(), highPrioClustersDR.end(), [&](auto& cluster) {
          communicationManager->progression();
          cluster->act();
        });
      }
      std::for_each(lowPrioClusters.begin(), lowPrioClusters.end(), [&](auto& cluster) {
        communicationManager->progression();
        cluster->act();
      });
      std::for_each(lowPrioClustersDR.begin(), lowPrioClustersDR.end(), [&](auto& cluster) {
        communicationManager->progression();
        cluster->act();
      });
      finished =
          std::all_of(clusters.begin(), clusters.end(), [](auto& c) { return c->synchronized(); });
      finished &= std::all_of(
          clustersDR.begin(), clustersDR.end(), [](auto& c) { return c->synchronized(); });
      finished &= communicationManager->checkIfFinished();
    }
  };

  std::invoke(phaseExecute);

#ifdef ACL_DEVICE
  device.api->syncDevice();
  device.api->popLastProfilingMark();
#endif

  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());
}

void TimeManager::printComputationTime(const std::string& outputPrefix,
                                       bool isLoopStatisticsNetcdfOutputOn) {
  actorStateStatisticsManager.finish();
  loopStatistics.printSummary(MPI::mpi.comm());
  loopStatistics.writeSamples(outputPrefix, isLoopStatisticsNetcdfOutputOn);
}

double TimeManager::getTimeTolerance() { return 1E-5 * timeStepping.globalCflTimeStepWidths[0]; }

void TimeManager::setPointSourcesForClusters(
    std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>>
        sourceClusters) {
  for (auto& cluster : clusters) {
    auto layerClusters = sourceClusters.find(cluster->getLayerType());
    if (layerClusters != sourceClusters.end() &&
        cluster->getClusterId() < layerClusters->second.size()) {
      cluster->setPointSources(std::move(layerClusters->second[cluster->getClusterId()]));
    }
  }
}

void TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter) {
  for (auto& cluster : clusters) {
    cluster->setReceiverCluster(
        receiverWriter.receiverCluster(cluster->getClusterId(), cluster->getLayerType()));
  }
}

void TimeManager::setInitialTimes(double time) {
  assert(time >= 0);

  for (auto& cluster : clusters) {
    cluster->setTime(time);
    cluster->setReceiverTime(time);
  }
}

void TimeManager::setTv(double tv) {
  for (auto& cluster : clusters) {
    cluster->setTv(tv);
  }
}

void TimeManager::freeDynamicResources() {
  for (auto& cluster : clusters) {
    cluster->freePointSources();
    cluster->finalize();
  }
  for (auto& cluster : clustersDR) {
    cluster->finalize();
  }
  communicationManager.reset(nullptr);
}

void TimeManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  const Executor exec = clusters[0]->getExecutor();
  bool sameExecutor = true;
  for (auto& cluster : clusters) {
    sameExecutor &= exec == cluster->getExecutor();
  }
  if (sameExecutor) {
    seissolInstance.getMemoryManager().synchronizeTo(place);
  } else {
    auto* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
    for (auto& cluster : clusters) {
      cluster->synchronizeTo(place, stream);
    }
    for (auto& cluster : clustersDR) {
      cluster->synchronizeTo(place, stream);
    }
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif
}

} // namespace seissol::time_stepping

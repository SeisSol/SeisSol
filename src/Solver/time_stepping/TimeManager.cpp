// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "Parallel/MPI.h"

#include "TimeManager.h"
#include "CommunicationManager.h"
#include "Monitoring/Instrumentation.h"
#include "SeisSol.h"
#include "ResultWriter/ClusteringWriter.h"
#include "Parallel/Helper.h"
#include <Memory/Tree/Layer.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif

seissol::time_stepping::TimeManager::TimeManager(seissol::SeisSol& seissolInstance):
  m_logUpdates(std::numeric_limits<unsigned int>::max()), seissolInstance(seissolInstance),
   actorStateStatisticsManager(m_loopStatistics)
{
  m_loopStatistics.addRegion("computeLocalIntegration");
  m_loopStatistics.addRegion("computeNeighboringIntegration");
  m_loopStatistics.addRegion("computeDynamicRupture");
  m_loopStatistics.addRegion("computePointSources");

  m_loopStatistics.enableSampleOutput(seissolInstance.getSeisSolParameters().output.loopStatisticsNetcdfOutput);
}

seissol::time_stepping::TimeManager::~TimeManager() {}

void seissol::time_stepping::TimeManager::addClusters(TimeStepping& timeStepping,
                                                      MeshStructure* i_meshStructure,
                                                      initializer::MemoryManager& memoryManager,
                                                      bool usePlasticity) {
  SCOREP_USER_REGION( "addClusters", SCOREP_USER_REGION_TYPE_FUNCTION );
  std::vector<std::unique_ptr<AbstractGhostTimeCluster>> ghostClusters;
  // assert non-zero pointers
  assert( i_meshStructure         != NULL );

  // store the time stepping
  m_timeStepping = timeStepping;

  auto clusteringWriter = writer::ClusteringWriter(memoryManager.getOutputPrefix());

  bool foundDynamicRuptureCluster = false;

  // iterate over local time clusters
  for (unsigned int localClusterId = 0; localClusterId < m_timeStepping.numberOfLocalClusters; localClusterId++) {
    // get memory layout of this cluster
    auto [meshStructure, globalData] = memoryManager.getMemoryLayout(localClusterId);

    const unsigned int l_globalClusterId = m_timeStepping.clusterIds[localClusterId];
    // chop off at synchronization time
    const auto timeStepSize = m_timeStepping.globalCflTimeStepWidths[l_globalClusterId];
    const long timeStepRate = ipow(static_cast<long>(m_timeStepping.globalTimeStepRates[0]),
         static_cast<long>(l_globalClusterId));

    // Dynamic rupture
    const auto interiorId = initializer::LayerIdentifier(HaloType::Interior, Config(), l_globalClusterId);
    const auto copyId = initializer::LayerIdentifier(HaloType::Copy, Config(), l_globalClusterId);
    const auto ghostId = initializer::LayerIdentifier(HaloType::Ghost, Config(), l_globalClusterId);

    // Note: We need to include the Ghost part, as we need to compute its DR part as well.
    const long numberOfDynRupCells
      = memoryManager.getDynamicRuptureTree()->layer(interiorId).size() +
        memoryManager.getDynamicRuptureTree()->layer(copyId).size() +
        memoryManager.getDynamicRuptureTree()->layer(ghostId).size();

    bool isFirstDynamicRuptureCluster = false;
    if (!foundDynamicRuptureCluster && numberOfDynRupCells > 0) {
      foundDynamicRuptureCluster = true;
      isFirstDynamicRuptureCluster = true;
    }
    auto& drScheduler = dynamicRuptureSchedulers.emplace_back(std::make_unique<DynamicRuptureScheduler>(numberOfDynRupCells,
                                                                                                        isFirstDynamicRuptureCluster));

    for (auto type : {HaloType::Copy, HaloType::Interior}) {
      const auto offsetMonitoring = type == HaloType::Interior ? 0 : m_timeStepping.numberOfGlobalClusters;
      // We print progress only if it is the cluster with the largest time step on each rank.
      // This does not mean that it is the largest cluster globally!
      const bool printProgress = (localClusterId == m_timeStepping.numberOfLocalClusters - 1) && (type == HaloType::Interior);
      const auto profilingId = l_globalClusterId + offsetMonitoring;
      auto* layerData = &memoryManager.getLtsTree()->layer(type == HaloType::Copy ? copyId : interiorId);
      auto* dynRupInteriorData = &memoryManager.getDynamicRuptureTree()->layer(interiorId);
      auto* dynRupCopyData = &memoryManager.getDynamicRuptureTree()->layer(copyId);
      clusters.push_back(std::make_unique<TimeCluster>(
          localClusterId,
          l_globalClusterId,
          profilingId,
          usePlasticity,
          layerData->getIdentifier(),
          timeStepSize,
          timeStepRate,
          printProgress,
          drScheduler.get(),
          globalData,
          layerData,
          dynRupInteriorData,
          dynRupCopyData,
          memoryManager.getLts(),
          memoryManager.getDynamicRupture(),
          memoryManager.getFrictionLaw(),
          memoryManager.getFrictionLawDevice(),
          memoryManager.getFaultOutputManager(),
          seissolInstance,
          &m_loopStatistics,
          &actorStateStatisticsManager.addCluster(profilingId))
      );

      const auto clusterSize = layerData->size();
      const auto dynRupSize = type == HaloType::Copy ? dynRupCopyData->size()
                                           : dynRupInteriorData->size();
      // Add writer to output
      clusteringWriter.addCluster(profilingId, localClusterId, type, clusterSize, dynRupSize);
    }
    auto& interior = clusters[clusters.size() - 1];
    auto& copy = clusters[clusters.size() - 2];

    // Mark copy layers as higher priority layers.
    interior->setPriority(ActorPriority::Low);
    copy->setPriority(ActorPriority::High);

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

    auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
    auto& [increaseManager, decreaseManager] = timeMirrorManagers;

    increaseManager.setTimeClusterVector(&clusters);
    decreaseManager.setTimeClusterVector(&clusters);

    // Create ghost time clusters for MPI
    const auto preferredDataTransferMode = MPI::mpi.getPreferredDataTransferMode();
    const auto persistent = usePersistentMpi(seissolInstance.env());
    const int globalClusterId = static_cast<int>(m_timeStepping.clusterIds[localClusterId]);
    for (unsigned int otherGlobalClusterId = 0; otherGlobalClusterId < m_timeStepping.numberOfGlobalClusters; ++otherGlobalClusterId) {
      const bool hasNeighborRegions = std::any_of(meshStructure->neighboringClusters,
                                                  meshStructure->neighboringClusters + meshStructure->numberOfRegions,
                                                  [otherGlobalClusterId](const auto& neighbor) {
        return static_cast<unsigned>(neighbor[1]) == otherGlobalClusterId;
      });
      if (hasNeighborRegions) {
        assert(static_cast<int>(otherGlobalClusterId) >= std::max(globalClusterId - 1, 0));
        assert(
            static_cast<int>(otherGlobalClusterId) <
            std::min(globalClusterId + 2, static_cast<int>(m_timeStepping.numberOfGlobalClusters)));
        const auto otherTimeStepSize = m_timeStepping.globalCflTimeStepWidths[otherGlobalClusterId];
        const long otherTimeStepRate =
            ipow(static_cast<long>(m_timeStepping.globalTimeStepRates[0]),
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
  }

  clusteringWriter.write();

  // Sort clusters by time step size in increasing order
  auto rateSorter = [](const auto& a, const auto& b) {
    return a->getTimeStepRate() < b->getTimeStepRate();
  };
  std::sort(clusters.begin(), clusters.end(), rateSorter);

  for (const auto& cluster : clusters) {
    if (cluster->getPriority() == ActorPriority::High) {
      highPrioClusters.emplace_back(cluster.get());
    } else {
      lowPrioClusters.emplace_back(cluster.get());
    }
  }

  std::sort(ghostClusters.begin(), ghostClusters.end(), rateSorter);

  if (seissol::useCommThread(MPI::mpi, seissolInstance.env())) {
    communicationManager = std::make_unique<ThreadedCommunicationManager>(std::move(ghostClusters),
                                                                          &seissolInstance.getPinning()
                                                                          );
  } else {
    communicationManager = std::make_unique<SerialCommunicationManager>(std::move(ghostClusters));
  }

  auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
  auto& [increaseManager, decreaseManager] = timeMirrorManagers;

  auto ghostClusterPointer = communicationManager->getGhostClusters();

  increaseManager.setGhostClusterVector(ghostClusterPointer);
  decreaseManager.setGhostClusterVector(ghostClusterPointer);

}

void seissol::time_stepping::TimeManager::setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
  m_faultOutputManager = faultOutputManager;
  for(auto& cluster : clusters) {
    cluster->setFaultOutputManager(faultOutputManager);
  }
}

seissol::dr::output::OutputManager* seissol::time_stepping::TimeManager::getFaultOutputManager() {
  assert(m_faultOutputManager != nullptr);
  return m_faultOutputManager;
}

void seissol::time_stepping::TimeManager::advanceInTime(const double &synchronizationTime) {
  SCOREP_USER_REGION( "advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION )

  // We should always move forward in time
  assert(m_timeStepping.synchronizationTime <= synchronizationTime);

  m_timeStepping.synchronizationTime = synchronizationTime;

  for (auto& cluster : clusters) {
    cluster->setSyncTime(synchronizationTime);
    cluster->reset();
  }

  communicationManager->reset(synchronizationTime);

  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());
#ifdef ACL_DEVICE
  device::DeviceInstance &device = device::DeviceInstance::getInstance();
  device.api->putProfilingMark("advanceInTime", device::ProfilingColors::Blue);
#endif

  // Move all clusters from RestartAfterSync to Corrected
  // Does not involve any computations
  for (auto& cluster : clusters) {
    assert(cluster->getNextLegalAction() == ActorAction::RestartAfterSync);
    cluster->act();
    assert(cluster->getState() == ActorState::Corrected);
  }

  bool finished = false; // Is true, once all clusters reached next sync point
  while (!finished) {
    finished = true;
    communicationManager->progression();

    // Update all high priority clusters
    std::for_each(highPrioClusters.begin(), highPrioClusters.end(), [&](auto& cluster) {
      if (cluster->getNextLegalAction() == ActorAction::Predict) {
        communicationManager->progression();
        cluster->act();
      }
    });
    std::for_each(highPrioClusters.begin(), highPrioClusters.end(), [&](auto& cluster) {
      if (cluster->getNextLegalAction() != ActorAction::Predict && cluster->getNextLegalAction() != ActorAction::Nothing) {
        communicationManager->progression();
        cluster->act();
      }
    });

    // Update one low priority cluster
    if (auto predictable = std::find_if(
          lowPrioClusters.begin(), lowPrioClusters.end(), [](auto& c) {
            return c->getNextLegalAction() == ActorAction::Predict;
          }
      );
        predictable != lowPrioClusters.end()) {
      (*predictable)->act();
    } else {
    }
    if (auto correctable = std::find_if(
          lowPrioClusters.begin(), lowPrioClusters.end(), [](auto& c) {
            return c->getNextLegalAction() != ActorAction::Predict && c->getNextLegalAction() != ActorAction::Nothing;
          }
      );
        correctable != lowPrioClusters.end()) {
      (*correctable)->act();
    } else {
    }
    finished = std::all_of(clusters.begin(), clusters.end(),
                           [](auto& c) {
      return c->synced();
    });
    finished &= communicationManager->checkIfFinished();
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
  for (auto& cluster : clusters) {
    cluster->finishPhase();
  }
}

void seissol::time_stepping::TimeManager::printComputationTime(
    const std::string& outputPrefix, bool isLoopStatisticsNetcdfOutputOn) {
  actorStateStatisticsManager.finish();
  m_loopStatistics.printSummary(MPI::mpi.comm());
  m_loopStatistics.writeSamples(outputPrefix, isLoopStatisticsNetcdfOutputOn);
}

double seissol::time_stepping::TimeManager::getTimeTolerance() {
  return 1E-5 * m_timeStepping.globalCflTimeStepWidths[0];
}

void seissol::time_stepping::TimeManager::setPointSourcesForClusters(
    std::vector<seissol::kernels::PointSourceClusterPair> sourceClusters) {
  const auto& map = seissolInstance.getMemoryManager().memoryContainer().colorMap;
  for (auto& cluster : clusters) {
    cluster->setPointSources(std::move(sourceClusters[map.colorId(cluster->getIdentifier())]));
  }
}

void seissol::time_stepping::TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter)
{
  for (auto& cluster : clusters) {
    cluster->setReceiverCluster(receiverWriter.receiverCluster(cluster->getIdentifier()));
  }
}

void seissol::time_stepping::TimeManager::setInitialTimes( double time ) {
  assert( time >= 0 );

  for (auto& cluster : clusters) {
    cluster->setPredictionTime(time);
    cluster->setCorrectionTime(time);
    cluster->setReceiverTime(time);
    cluster->setLastSubTime(time);
  }
  for (auto& cluster : *communicationManager->getGhostClusters()) {
    cluster->setPredictionTime(time);
    cluster->setCorrectionTime(time);
  }
}

void seissol::time_stepping::TimeManager::setTv(double tv) {
  for(auto& cluster : clusters) {
    cluster->setTv(tv);
  }
}

void seissol::time_stepping::TimeManager::freeDynamicResources() {
  for (auto& cluster : clusters) {
    cluster->freePointSources();
    cluster->finalize();
  }
  communicationManager.reset(nullptr);
}

void seissol::time_stepping::TimeManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  Executor exec = clusters[0]->getExecutor();
  bool sameExecutor = true;
  for (auto& cluster : clusters) {
    sameExecutor &= exec == cluster->getExecutor();
  }
  if (sameExecutor) {
    seissolInstance.getMemoryManager().synchronizeTo(place);
  }
  else {
    auto* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
    for (auto& cluster : clusters) {
      cluster->synchronizeTo(place, stream);
    }
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif
}


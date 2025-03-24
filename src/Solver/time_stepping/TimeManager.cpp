// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
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
#include "Initializer/TimeStepping/Common.h"
#include "SeisSol.h"
#include "ResultWriter/ClusteringWriter.h"
#include "Parallel/Helper.h"
#include "Solver/MultipleSimulations.h"

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

seissol::time_stepping::TimeManager::~TimeManager() = default;

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

    std::array<long, seissol::multipleSimulations::numberOfSimulations> numDynRupCells;
    // Dynamic rupture
    for (unsigned i = 0; i < seissol::multipleSimulations::numberOfSimulations; i++) {
      auto& dynRupTree = memoryManager.getDynamicRuptureTree()[i]->child(localClusterId);
      const long numberOfDynRupCells = dynRupTree.child(Interior).getNumberOfCells() +
        dynRupTree.child(Copy).getNumberOfCells() +
        dynRupTree.child(Ghost).getNumberOfCells();
      numDynRupCells[i] = numberOfDynRupCells;
    } // (TO DISCUSS: does this make remotely any sense? Do all the LTS Trees for Dynamic Rupture have the same number of Cells in Copy, Ghost etc.? If so,
    // can I do this just once and not run it for the entire loop? I don't understand why we need an array if that is the case)

    // Note: We need to include the Ghost part, as we need to compute its DR part as well.
 
    bool isFirstDynamicRuptureCluster = false;
    if (!foundDynamicRuptureCluster && numDynRupCells[0] > 0) {
      foundDynamicRuptureCluster = true;
      isFirstDynamicRuptureCluster = true;
    }

    auto& drScheduler = dynamicRuptureSchedulers.emplace_back(std::make_unique<DynamicRuptureScheduler>(numDynRupCells[0],
                                                                                                        isFirstDynamicRuptureCluster));
                                                                                                        //(TO DISCUSS: Same as the above issue with the dynRupTree)

    for (auto type : {Copy, Interior}) {
      const auto offsetMonitoring = type == Interior ? 0 : m_timeStepping.numberOfGlobalClusters;
      // We print progress only if it is the cluster with the largest time step on each rank.
      // This does not mean that it is the largest cluster globally!
      const bool printProgress = (localClusterId == m_timeStepping.numberOfLocalClusters - 1) && (type == Interior);
      const auto profilingId = l_globalClusterId + offsetMonitoring;
      auto* layerData = &memoryManager.getLtsTree()->child(localClusterId).child(type);
      std::array<seissol::initializer::Layer*, seissol::multipleSimulations::numberOfSimulations> dynRupInteriorData;
      std::array<seissol::initializer::Layer*, seissol::multipleSimulations::numberOfSimulations> dynRupCopyData;

      // auto* dynRupInteriorData = &dynRupTree.child(Interior);
      // auto* dynRupCopyData = &dynRupTree.child(Copy);
      for (unsigned int i = 0; i < seissol::multipleSimulations::numberOfSimulations; i++) {
        dynRupCopyData[i] =
            &memoryManager.getDynamicRuptureTree()[i]->child(localClusterId).child(Copy);
        dynRupInteriorData[i] =
            &memoryManager.getDynamicRuptureTree()[i]->child(localClusterId).child(Interior);
      } // (TO DISCUSS: Do we need different LTS Trees just because we would need different Copy
        // datas and Interior Datas? Can't we change the structure of the LTS Tree to get multiple
        // of these there?)

      clusters.push_back(std::make_unique<TimeCluster>(
          localClusterId,
          l_globalClusterId,
          profilingId,
          usePlasticity,
          type,
          timeStepSize,
          timeStepRate,
          printProgress,
          drScheduler.get(), // ignore it
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

      const auto clusterSize = layerData->getNumberOfCells();
      const auto dynRupSize =
          type == Copy ? dynRupCopyData[0]->getNumberOfCells()
                       : dynRupInteriorData[0]
                             ->getNumberOfCells(); // (TO DISCUSS: Confirmation that number of cells
                                                   // would remain the same and that miniseissol
                                                   // won't do some magic to give different number
                                                   // of cells for different simulations)
      // Add writer to output
      clusteringWriter.addCluster(profilingId, localClusterId, type, clusterSize, dynRupSize);
      //ignore it for now
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
#ifdef USE_MPI
    // Create ghost time clusters for MPI
    const auto preferredDataTransferMode = MPI::mpi.getPreferredDataTransferMode();
    const auto persistent = usePersistentMpi();
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
#endif
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

  if (seissol::useCommThread(MPI::mpi)) {
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

void seissol::time_stepping::TimeManager::setFaultOutputManager(
    std::array<std::shared_ptr<seissol::dr::output::OutputManager>,
               seissol::multipleSimulations::numberOfSimulations>& faultOutputManager) {
  for (unsigned int i = 0; i < seissol::multipleSimulations::numberOfSimulations; i++) {
    m_faultOutputManager[i] = faultOutputManager[i].get();
  }

  for (auto& cluster : clusters) {
    cluster->setFaultOutputManager(faultOutputManager);
  }
}

std::array<seissol::dr::output::OutputManager*, seissol::multipleSimulations::numberOfSimulations>&
    seissol::time_stepping::TimeManager::getFaultOutputManager() {
  for (unsigned int i = 0; i < seissol::multipleSimulations::numberOfSimulations; i++) {
    assert(m_faultOutputManager[i] != nullptr);
  }
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
    std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>> sourceClusters) {
  for (auto& cluster : clusters) {
    auto layerClusters = sourceClusters.find(cluster->getLayerType());
    if (layerClusters != sourceClusters.end() && cluster->getClusterId() < layerClusters->second.size()) {
      cluster->setPointSources(std::move(layerClusters->second[cluster->getClusterId()]));
    }
  }
}

void seissol::time_stepping::TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter)
{
  for (auto& cluster : clusters) {
    cluster->setReceiverCluster(receiverWriter.receiverCluster(cluster->getClusterId(),
                                                               cluster->getLayerType()));
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


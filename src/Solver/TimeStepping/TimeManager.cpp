// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "Parallel/MPI.h"

#include "CommunicationManager.h"
#include "Parallel/Helper.h"
#include "ResultWriter/ClusteringWriter.h"
#include "SeisSol.h"
#include "TimeManager.h"
#include <DynamicRupture/Output/OutputManager.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Typedefs.h>
#include <Kernels/PointSourceCluster.h>
#include <Memory/Tree/Layer.h>
#include <Monitoring/Instrumentation.h>
#include <ResultWriter/ReceiverWriter.h>
#include <Solver/TimeStepping/AbstractGhostTimeCluster.h>
#include <Solver/TimeStepping/AbstractTimeCluster.h>
#include <Solver/TimeStepping/ActorState.h>
#include <Solver/TimeStepping/GhostTimeClusterFactory.h>
#include <Solver/TimeStepping/TimeCluster.h>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::time_stepping {

TimeManager::TimeManager(seissol::SeisSol& seissolInstance)
    : seissolInstance(seissolInstance), actorStateStatisticsManager(loopStatistics) {
  loopStatistics.addRegion("computeLocalIntegration");
  loopStatistics.addRegion("computeNeighboringIntegration");
  loopStatistics.addRegion("computeDynamicRupture");
  loopStatistics.addRegion("computePointSources");

  loopStatistics.enableSampleOutput(
      seissolInstance.getSeisSolParameters().output.loopStatisticsNetcdfOutput);
}

TimeManager::~TimeManager() = default;

void TimeManager::addClusters(const initializer::ClusterLayout& clusterLayout,
                              MeshStructure* meshStructure,
                              initializer::MemoryManager& memoryManager,
                              bool usePlasticity) {
  SCOREP_USER_REGION("addClusters", SCOREP_USER_REGION_TYPE_FUNCTION);
  std::vector<std::unique_ptr<AbstractGhostTimeCluster>> ghostClusters;
  // assert non-zero pointers
  assert(meshStructure != nullptr);

  // store the time stepping
  this->clusterLayout = clusterLayout;

  auto clusteringWriter = writer::ClusteringWriter(memoryManager.getOutputPrefix());

  bool foundDynamicRuptureCluster = false;

  // iterate over local time clusters
  for (std::size_t localClusterId = 0; localClusterId < clusterLayout.globalClusterCount;
       ++localClusterId) {
    // get memory layout of this cluster
    auto [meshStructure, globalData] = memoryManager.getMemoryLayout(localClusterId);

    const int globalClusterId = static_cast<int>(localClusterId);
    // chop off at synchronization time
    const auto timeStepSize = clusterLayout.timestepRate(localClusterId);
    const auto timeStepRate = clusterLayout.clusterRate(localClusterId);

    // Dynamic rupture
    auto& dynRupTree = memoryManager.getDynamicRuptureTree()->child(localClusterId);
    // Note: We need to include the Ghost part, as we need to compute its DR part as well.
    const long numberOfDynRupCells = dynRupTree.child(Interior).size() +
                                     dynRupTree.child(Copy).size() + dynRupTree.child(Ghost).size();

    bool isFirstDynamicRuptureCluster = false;
    if (!foundDynamicRuptureCluster && numberOfDynRupCells > 0) {
      foundDynamicRuptureCluster = true;
      isFirstDynamicRuptureCluster = true;
    }
    auto& drScheduler =
        dynamicRuptureSchedulers.emplace_back(std::make_unique<DynamicRuptureScheduler>(
            numberOfDynRupCells, isFirstDynamicRuptureCluster));

    for (auto type : {Copy, Interior}) {
      const auto offsetMonitoring = type == Interior ? 0 : clusterLayout.globalClusterCount;
      // We print progress only if it is the cluster with the largest time step on each rank.
      // This does not mean that it is the largest cluster globally!
      const bool printProgress =
          (localClusterId == clusterLayout.globalClusterCount - 1) && (type == Interior);
      const auto profilingId = globalClusterId + offsetMonitoring;
      auto* layerData = &memoryManager.getLtsTree()->child(localClusterId).child(type);
      auto* dynRupInteriorData = &dynRupTree.child(Interior);
      auto* dynRupCopyData = &dynRupTree.child(Copy);
      clusters.push_back(
          std::make_unique<TimeCluster>(localClusterId,
                                        globalClusterId,
                                        profilingId,
                                        usePlasticity,
                                        type,
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
                                        &loopStatistics,
                                        &actorStateStatisticsManager.addCluster(profilingId)));

      const auto clusterSize = layerData->size();
      const auto dynRupSize = type == Copy ? dynRupCopyData->size() : dynRupInteriorData->size();
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
      for (int i = 0; i < 2; ++i) {
        copy->connect(*clusters[clusters.size() - 2 - i - 1]);
        interior->connect(*clusters[clusters.size() - 2 - i - 1]);
      }
    }

    // Create ghost time clusters for MPI
    const auto preferredDataTransferMode = MPI::mpi.getPreferredDataTransferMode();
    const auto persistent = usePersistentMpi(seissolInstance.env());
    for (unsigned int otherGlobalClusterId = 0;
         otherGlobalClusterId < clusterLayout.globalClusterCount;
         ++otherGlobalClusterId) {
      const bool hasNeighborRegions =
          std::any_of(meshStructure->neighboringClusters,
                      meshStructure->neighboringClusters + meshStructure->numberOfRegions,
                      [otherGlobalClusterId](const auto& neighbor) {
                        return static_cast<unsigned>(neighbor[1]) == otherGlobalClusterId;
                      });
      if (hasNeighborRegions) {
        assert(static_cast<int>(otherGlobalClusterId) >= std::max(globalClusterId - 1, 0));
        assert(static_cast<int>(otherGlobalClusterId) <
               std::min(globalClusterId + 2, static_cast<int>(clusterLayout.globalClusterCount)));
        const auto otherTimeStepSize = clusterLayout.timestepRate(otherGlobalClusterId);
        const auto otherTimeStepRate = clusterLayout.clusterRate(otherGlobalClusterId);

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
    communicationManager = std::make_unique<ThreadedCommunicationManager>(
        std::move(ghostClusters), &seissolInstance.getPinning());
  } else {
    communicationManager = std::make_unique<SerialCommunicationManager>(std::move(ghostClusters));
  }

  auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
  auto& [increaseManager, decreaseManager] = timeMirrorManagers;

  auto* ghostClusterPointer = communicationManager->getGhostClusters();

  std::vector<AbstractTimeCluster*> allClusters(clusters.size() + ghostClusterPointer->size());
  for (std::size_t i = 0; i < clusters.size(); ++i) {
    allClusters[i] = clusters[i].get();
  }
  for (std::size_t i = 0; i < ghostClusterPointer->size(); ++i) {
    allClusters[i + clusters.size()] = ghostClusterPointer->at(i).get();
  }

  increaseManager.setClusterVector(allClusters);
  decreaseManager.setClusterVector(allClusters);
}

void TimeManager::setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
  this->faultOutputManager = faultOutputManager;
  for (auto& cluster : clusters) {
    cluster->setFaultOutputManager(faultOutputManager);
  }
}

seissol::dr::output::OutputManager* TimeManager::getFaultOutputManager() {
  assert(faultOutputManager != nullptr);
  return faultOutputManager;
}

void TimeManager::advanceInTime(const double& synchronizationTime) {
  SCOREP_USER_REGION("advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION)

  for (auto& cluster : clusters) {
    cluster->setSyncTime(synchronizationTime);
    cluster->reset();
  }

  communicationManager->reset(synchronizationTime);

  seissol::MPI::barrier(seissol::MPI::mpi.comm());
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
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
    communicationManager->progression();

    // Update all high priority clusters
    std::for_each(highPrioClusters.begin(), highPrioClusters.end(), [&](auto& cluster) {
      if (cluster->getNextLegalAction() == ActorAction::Predict) {
        communicationManager->progression();
        cluster->act();
      }
    });
    std::for_each(highPrioClusters.begin(), highPrioClusters.end(), [&](auto& cluster) {
      if (cluster->getNextLegalAction() != ActorAction::Predict &&
          cluster->getNextLegalAction() != ActorAction::Nothing) {
        communicationManager->progression();
        cluster->act();
      }
    });

    // Update one low priority cluster
    if (auto predictable =
            std::find_if(lowPrioClusters.begin(),
                         lowPrioClusters.end(),
                         [](auto& c) { return c->getNextLegalAction() == ActorAction::Predict; });
        predictable != lowPrioClusters.end()) {
      (*predictable)->act();
    } else {
    }
    if (auto correctable = std::find_if(lowPrioClusters.begin(),
                                        lowPrioClusters.end(),
                                        [](auto& c) {
                                          return c->getNextLegalAction() != ActorAction::Predict &&
                                                 c->getNextLegalAction() != ActorAction::Nothing;
                                        });
        correctable != lowPrioClusters.end()) {
      (*correctable)->act();
    } else {
    }
    finished = std::all_of(clusters.begin(), clusters.end(), [](auto& c) { return c->synced(); });
    finished &= communicationManager->checkIfFinished();
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
  for (auto& cluster : clusters) {
    cluster->finishPhase();
  }
}

void TimeManager::printComputationTime(const std::string& outputPrefix,
                                       bool isLoopStatisticsNetcdfOutputOn) {
  actorStateStatisticsManager.finish();
  loopStatistics.printSummary(MPI::mpi.comm());
  loopStatistics.writeSamples(outputPrefix, isLoopStatisticsNetcdfOutputOn);
}

double TimeManager::getTimeTolerance() const {
  return 1E-5 * clusterLayout.value().minimumTimestep;
}

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
  }
  for (auto& cluster : *communicationManager->getGhostClusters()) {
    cluster->setTime(time);
  }
}

void TimeManager::freeDynamicResources() {
  for (auto& cluster : clusters) {
    cluster->finalize();
  }
  communicationManager.reset(nullptr);
}

void TimeManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  Executor exec = clusters[0]->getExecutor();
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
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif
}

} // namespace seissol::time_stepping

// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "TimeManager.h"

#include "Common/Iterator.h"
#include "CommunicationManager.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Kernels/PointSourceCluster.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "ResultWriter/ClusteringWriter.h"
#include "ResultWriter/ReceiverWriter.h"
#include "SeisSol.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/GhostTimeClusterFactory.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/TimeCluster.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <mpi.h>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
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
                              const solver::HaloCommunication& haloStructure,
                              initializer::MemoryManager& memoryManager,
                              bool usePlasticity) {
  SCOREP_USER_REGION("addClusters", SCOREP_USER_REGION_TYPE_FUNCTION);
  std::vector<std::unique_ptr<AbstractGhostTimeCluster>> ghostClusters;

  // store the time stepping
  this->clusterLayout = clusterLayout;

  auto clusteringWriter =
      writer::ClusteringWriter(seissolInstance.getSeisSolParameters().output.prefix);

  std::vector<std::size_t> drCellsPerCluster(clusterLayout.globalClusterCount);

  // setup DR schedulers
  for (const auto& layer : memoryManager.getDRStorage().leaves()) {
    drCellsPerCluster[layer.getIdentifier().lts] += layer.size();
  }

  std::size_t drClusterOutput = std::numeric_limits<std::size_t>::max();
  for (std::size_t clusterId = 0; clusterId < drCellsPerCluster.size(); ++clusterId) {
    if (drCellsPerCluster[clusterId] > 0) {
      drClusterOutput = clusterId;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE,
                &drClusterOutput,
                1,
                MPI::castToMpiType<std::size_t>(),
                MPI_MIN,
                MPI::mpi.comm());

  for (std::size_t clusterId = 0; clusterId < drCellsPerCluster.size(); ++clusterId) {
    const bool isFirstDynamicRuptureCluster = drClusterOutput == clusterId;
    dynamicRuptureSchedulers.emplace_back(std::make_unique<DynamicRuptureScheduler>(
        drCellsPerCluster[clusterId], isFirstDynamicRuptureCluster));
  }

  std::vector<AbstractTimeCluster*> cellClusterBackmap(
      memoryManager.getLtsStorage().getColorMap().size());

  const auto deltaId = [&](const auto& id, HaloType halo, int32_t offset) {
    auto cloned = id;
    cloned.halo = halo;
    cloned.lts += offset;
    return memoryManager.getLtsStorage().getColorMap().colorId(cloned);
  };

  // iterate over local time clusters
  for (auto& layer : memoryManager.getLtsStorage().leaves(Ghost)) {
    auto globalData = memoryManager.getGlobalData();

    const auto clusterId = layer.getIdentifier().lts;

    // chop off at synchronization time
    const auto timeStepSize = clusterLayout.timestepRate(clusterId);
    const auto timeStepRate = clusterLayout.clusterRate(clusterId);

    // We print progress only if it is the cluster with the largest time step on each rank.
    // This does not mean that it is the largest cluster globally!
    const bool printProgress = (clusterId == clusterLayout.globalClusterCount - 1) &&
                               (layer.getIdentifier().halo == HaloType::Interior);
    const auto profilingId = layer.id();

    auto* dynRupInteriorData =
        &memoryManager.getDRStorage().layer(deltaId(layer.getIdentifier(), HaloType::Interior, 0));
    auto* dynRupCopyData =
        &memoryManager.getDRStorage().layer(deltaId(layer.getIdentifier(), HaloType::Copy, 0));

    auto& cluster = clusters.emplace_back(
        std::make_unique<TimeCluster>(clusterId,
                                      clusterId,
                                      profilingId,
                                      usePlasticity,
                                      layer.getIdentifier().halo,
                                      timeStepSize,
                                      timeStepRate,
                                      printProgress,
                                      dynamicRuptureSchedulers[clusterId].get(),
                                      globalData,
                                      &layer,
                                      dynRupInteriorData,
                                      dynRupCopyData,
                                      memoryManager.getFrictionLaw(),
                                      memoryManager.getFrictionLawDevice(),
                                      memoryManager.getFaultOutputManager(),
                                      seissolInstance,
                                      &loopStatistics,
                                      &actorStateStatisticsManager.addCluster(profilingId)));

    const auto clusterSize = layer.size();
    const auto dynRupSize = memoryManager.getDRStorage().layer(layer.id()).size();
    // Add writer to output
    clusteringWriter.addCluster(
        profilingId, clusterId, layer.getIdentifier().halo, clusterSize, dynRupSize);

    if (layer.getIdentifier().halo == HaloType::Copy) {
      cluster->setPriority(ActorPriority::High);
    } else {
      cluster->setPriority(ActorPriority::Low);
    }

    cellClusterBackmap[layer.id()] = cluster.get();
  }

  const auto connectIfBothExist = [](auto& a, auto& b) {
    if (a != nullptr && b != nullptr) {
      a->connect(*b);
    }
  };

  for (auto& layer : memoryManager.getLtsStorage().leaves(Ghost)) {
    for (auto& otherLayer : memoryManager.getLtsStorage().leaves(Ghost)) {
      // only traverse half of all combinations
      if (layer.id() < otherLayer.id()) {
        const int64_t lts1 = layer.getIdentifier().lts;
        const int64_t lts2 = otherLayer.getIdentifier().lts;
        if (std::abs(lts1 - lts2) <= 1) {
          connectIfBothExist(cellClusterBackmap[layer.id()], cellClusterBackmap[otherLayer.id()]);
        }
      }
    }
  }

  // Create ghost time clusters for MPI
  const auto preferredDataTransferMode = MPI::mpi.getPreferredDataTransferMode();
  const auto persistent = usePersistentMpi(seissolInstance.env());
  for (auto& layer : memoryManager.getLtsStorage().leaves(Ghost | Interior)) {

    const auto clusterId = layer.getIdentifier().lts;

    for (const auto [i, halo] : common::enumerate(haloStructure.at(layer.id()))) {

      const bool hasNeighborRegions = !halo.copy.empty() || !halo.ghost.empty();
      const auto other = memoryManager.getLtsStorage().getColorMap().argument(i);

      if (hasNeighborRegions) {

        assert(other.halo == HaloType::Ghost);
        assert(other.lts + 1 >= clusterId);
        assert(other.lts < std::min(clusterId + 2, clusterLayout.globalClusterCount));

        const auto otherTimeStepSize = clusterLayout.timestepRate(other.lts);
        const auto otherTimeStepRate = clusterLayout.clusterRate(other.lts);

        auto ghostCluster = GhostTimeClusterFactory::get(otherTimeStepSize,
                                                         otherTimeStepRate,
                                                         layer.id(),
                                                         i,
                                                         haloStructure,
                                                         preferredDataTransferMode,
                                                         persistent);
        ghostClusters.push_back(std::move(ghostCluster));

        // Connect with previous copy layer.
        ghostClusters.back()->connect(*cellClusterBackmap[layer.id()]);
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

    // dereference first due to clang-tidy recommendation
    (*cluster).reset();
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
    std::vector<seissol::kernels::PointSourceClusterPair> sourceClusters) {
  for (auto& cluster : clusters) {
    cluster->setPointSources(std::move(sourceClusters[cluster->layerId()]));
  }
}

void TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter) {
  for (auto& cluster : clusters) {
    cluster->setReceiverCluster(receiverWriter.receiverCluster(cluster->layerId()));
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

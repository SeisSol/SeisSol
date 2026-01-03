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
    : seissolInstance_(seissolInstance), actorStateStatisticsManager_(loopStatistics_) {
  loopStatistics_.addRegion("computeLocalIntegration");
  loopStatistics_.addRegion("computeNeighboringIntegration");
  loopStatistics_.addRegion("computeDynamicRupture");
  loopStatistics_.addRegion("computePointSources");

  loopStatistics_.enableSampleOutput(
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
  this->clusterLayout_ = clusterLayout;

  auto clusteringWriter =
      writer::ClusteringWriter(seissolInstance_.getSeisSolParameters().output.prefix);

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
                Mpi::castToMpiType<std::size_t>(),
                MPI_MIN,
                Mpi::mpi.comm());

  for (std::size_t clusterId = 0; clusterId < drCellsPerCluster.size(); ++clusterId) {
    const bool isFirstDynamicRuptureCluster = drClusterOutput == clusterId;
    dynamicRuptureSchedulers_.emplace_back(std::make_unique<DynamicRuptureScheduler>(
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

    auto& cluster = clusters_.emplace_back(
        std::make_unique<TimeCluster>(clusterId,
                                      clusterId,
                                      profilingId,
                                      usePlasticity,
                                      layer.getIdentifier().halo,
                                      timeStepSize,
                                      timeStepRate,
                                      printProgress,
                                      dynamicRuptureSchedulers_[clusterId].get(),
                                      globalData,
                                      &layer,
                                      dynRupInteriorData,
                                      dynRupCopyData,
                                      memoryManager.getFrictionLaw(),
                                      memoryManager.getFrictionLawDevice(),
                                      memoryManager.getFaultOutputManager(),
                                      seissolInstance_,
                                      &loopStatistics_,
                                      &actorStateStatisticsManager_.addCluster(profilingId)));

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
  const auto preferredDataTransferMode = Mpi::mpi.getPreferredDataTransferMode();
  const auto persistent = usePersistentMpi(seissolInstance_.env());
  for (auto& layer : memoryManager.getLtsStorage().leaves(Ghost | Interior)) {

    for (const auto [i, halo] : common::enumerate(haloStructure.at(layer.id()))) {

      const bool hasNeighborRegions = !halo.copy.empty() || !halo.ghost.empty();
      const auto other = memoryManager.getLtsStorage().getColorMap().argument(i);

      if (hasNeighborRegions) {

        assert(other.halo == HaloType::Ghost);
        assert(other.lts + 1 >= layer.getIdentifier().lts);
        assert(other.lts <
               std::min(layer.getIdentifier().lts + 2, clusterLayout.globalClusterCount));

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
  std::sort(clusters_.begin(), clusters_.end(), rateSorter);

  for (const auto& cluster : clusters_) {
    if (cluster->getPriority() == ActorPriority::High) {
      highPrioClusters_.emplace_back(cluster.get());
    } else {
      lowPrioClusters_.emplace_back(cluster.get());
    }
  }

  std::sort(ghostClusters.begin(), ghostClusters.end(), rateSorter);

  if (seissol::useCommThread(Mpi::mpi, seissolInstance_.env())) {
    communicationManager_ = std::make_unique<ThreadedCommunicationManager>(
        std::move(ghostClusters), &seissolInstance_.getPinning());
  } else {
    communicationManager_ = std::make_unique<SerialCommunicationManager>(std::move(ghostClusters));
  }

  auto& timeMirrorManagers = seissolInstance_.getTimeMirrorManagers();
  auto& [increaseManager, decreaseManager] = timeMirrorManagers;

  auto* ghostClusterPointer = communicationManager_->getGhostClusters();

  std::vector<AbstractTimeCluster*> allClusters(clusters_.size() + ghostClusterPointer->size());
  for (std::size_t i = 0; i < clusters_.size(); ++i) {
    allClusters[i] = clusters_[i].get();
  }
  for (std::size_t i = 0; i < ghostClusterPointer->size(); ++i) {
    allClusters[i + clusters_.size()] = ghostClusterPointer->at(i).get();
  }

  increaseManager.setClusterVector(allClusters);
  decreaseManager.setClusterVector(allClusters);
}

void TimeManager::setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
  this->faultOutputManager_ = faultOutputManager;
  for (auto& cluster : clusters_) {
    cluster->setFaultOutputManager(faultOutputManager);
  }
}

seissol::dr::output::OutputManager* TimeManager::getFaultOutputManager() {
  assert(faultOutputManager != nullptr);
  return faultOutputManager_;
}

void TimeManager::advanceInTime(const double& synchronizationTime) {
  SCOREP_USER_REGION("advanceInTime", SCOREP_USER_REGION_TYPE_FUNCTION)

  for (auto& cluster : clusters_) {
    cluster->setSyncTime(synchronizationTime);

    // dereference first due to clang-tidy recommendation
    (*cluster).reset();
  }

  communicationManager_->reset(synchronizationTime);

  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->putProfilingMark("advanceInTime", device::ProfilingColors::Blue);
#endif

  // Move all clusters from RestartAfterSync to Corrected
  // Does not involve any computations
  for (auto& cluster : clusters_) {
    assert(cluster->getNextLegalAction() == ActorAction::RestartAfterSync);
    cluster->act();
    assert(cluster->getState() == ActorState::Corrected);
  }

  bool finished = false; // Is true, once all clusters reached next sync point
  while (!finished) {
    communicationManager_->progression();

    // Update all high priority clusters
    std::for_each(highPrioClusters_.begin(), highPrioClusters_.end(), [&](auto& cluster) {
      if (cluster->getNextLegalAction() == ActorAction::Predict) {
        communicationManager_->progression();
        cluster->act();
      }
    });
    std::for_each(highPrioClusters_.begin(), highPrioClusters_.end(), [&](auto& cluster) {
      if (cluster->getNextLegalAction() != ActorAction::Predict &&
          cluster->getNextLegalAction() != ActorAction::Nothing) {
        communicationManager_->progression();
        cluster->act();
      }
    });

    // Update one low priority cluster
    if (auto predictable =
            std::find_if(lowPrioClusters_.begin(),
                         lowPrioClusters_.end(),
                         [](auto& c) { return c->getNextLegalAction() == ActorAction::Predict; });
        predictable != lowPrioClusters_.end()) {
      (*predictable)->act();
    } else {
    }
    if (auto correctable = std::find_if(lowPrioClusters_.begin(),
                                        lowPrioClusters_.end(),
                                        [](auto& c) {
                                          return c->getNextLegalAction() != ActorAction::Predict &&
                                                 c->getNextLegalAction() != ActorAction::Nothing;
                                        });
        correctable != lowPrioClusters_.end()) {
      (*correctable)->act();
    } else {
    }
    finished = std::all_of(clusters_.begin(), clusters_.end(), [](auto& c) { return c->synced(); });
    finished &= communicationManager_->checkIfFinished();
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
  for (auto& cluster : clusters_) {
    cluster->finishPhase();
  }
}

void TimeManager::printComputationTime(const std::string& outputPrefix,
                                       bool isLoopStatisticsNetcdfOutputOn) {
  actorStateStatisticsManager_.finish();
  loopStatistics_.printSummary(Mpi::mpi.comm());
  loopStatistics_.writeSamples(outputPrefix, isLoopStatisticsNetcdfOutputOn);
}

double TimeManager::getTimeTolerance() const {
  return 1E-5 * clusterLayout_.value().minimumTimestep;
}

void TimeManager::setPointSourcesForClusters(
    std::vector<seissol::kernels::PointSourceClusterPair> sourceClusters) {
  for (auto& cluster : clusters_) {
    cluster->setPointSources(std::move(sourceClusters[cluster->layerId()]));
  }
}

void TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter) {
  for (auto& cluster : clusters_) {
    cluster->setReceiverCluster(receiverWriter.receiverCluster(cluster->layerId()));
  }
}

void TimeManager::setInitialTimes(double time) {
  assert(time >= 0);

  for (auto& cluster : clusters_) {
    cluster->setTime(time);
  }
  for (auto& cluster : *communicationManager_->getGhostClusters()) {
    cluster->setTime(time);
  }
}

void TimeManager::freeDynamicResources() {
  for (auto& cluster : clusters_) {
    cluster->finalize();
  }
  communicationManager_.reset(nullptr);
}

void TimeManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  const auto exec = clusters_[0]->getExecutor();
  bool sameExecutor = true;
  for (auto& cluster : clusters_) {
    sameExecutor &= exec == cluster->getExecutor();
  }
  if (sameExecutor) {
    seissolInstance_.getMemoryManager().synchronizeTo(place);
  } else {
    auto* stream = device::DeviceInstance::getInstance().api->getDefaultStream();
    for (auto& cluster : clusters_) {
      cluster->synchronizeTo(place, stream);
    }
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif
}

} // namespace seissol::time_stepping

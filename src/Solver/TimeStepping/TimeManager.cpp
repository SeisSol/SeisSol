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
#include "Kernels/Common.h"
#include "Kernels/PointSourceCluster.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "ResultWriter/ClusteringWriter.h"
#include "ResultWriter/ReceiverWriter.h"
#include "SeisSol.h"
#include "Solver/Settings.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/CellCluster.h"
#include "Solver/TimeStepping/GhostCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/HaloTransport.h"
#include "Solver/TimeStepping/TimeSteppingPlan.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/AbstractAPI.h>
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
      seissolInstance.parameters().output.loopStatisticsNetcdfOutput);
}

TimeManager::~TimeManager() = default;

namespace {

/**
 * Compares the connections between the clusters with the dependencies that the face neighbors of
 * the cells impose: each cell layer is connected to every local layer of a time cluster differing
 * by at most one, and a copy layer additionally to one ghost cluster per ghost layer it exchanges
 * data with. Dependencies without a connection are an error; connections without a dependency
 * only cost synchronization.
 */
void reportClusterDependencies(initializer::MemoryManager& memoryManager,
                               const solver::HaloCommunication& haloStructure) {
  const auto& colorMap = memoryManager.ltsStorage().getColorMap();

  std::size_t connections = 0;
  std::size_t dependencies = 0;
  std::size_t missing = 0;
  for (auto& layer : memoryManager.ltsStorage().leaves(Ghost)) {
    std::vector<bool> connected(colorMap.size(), false);
    for (const auto& other : memoryManager.ltsStorage().leaves(Ghost)) {
      const auto lts1 = static_cast<int64_t>(layer.getIdentifier().lts);
      const auto lts2 = static_cast<int64_t>(other.getIdentifier().lts);
      connected[other.id()] = other.id() != layer.id() && std::abs(lts1 - lts2) <= 1;
    }
    if (layer.getIdentifier().halo == HaloType::Copy) {
      for (const auto [color, halo] : common::enumerate(haloStructure.at(layer.id()))) {
        connected[color] = connected[color] || !halo.copy.empty() || !halo.ghost.empty();
      }
    }

    std::vector<bool> needed(colorMap.size(), false);
    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (const auto& neighbor : secondaryInformation[cell].faceNeighbors) {
        if (neighbor.color < colorMap.size() && neighbor.color != layer.id()) {
          needed[neighbor.color] = true;
        }
      }
    }

    for (std::size_t color = 0; color < colorMap.size(); ++color) {
      connections += connected[color] ? 1 : 0;
      dependencies += needed[color] ? 1 : 0;
      missing += (needed[color] && !connected[color]) ? 1 : 0;
    }
  }

  std::array<std::size_t, 3> counts{connections, dependencies, missing};
  MPI_Allreduce(MPI_IN_PLACE,
                counts.data(),
                counts.size(),
                Mpi::castToMpiType<std::size_t>(),
                MPI_SUM,
                Mpi::mpi.comm());

  logInfo() << "Cell cluster connections:" << counts[0]
            << "; needed by face neighbors:" << counts[1] << "(summed over all ranks)";
  if (counts[2] > 0) {
    logError() << counts[2] << "dependencies between cell layers have no connection.";
  }
}

} // namespace

void TimeManager::addClusters(const initializer::ClusterLayout& clusterLayout,
                              const solver::HaloCommunication& haloStructure,
                              initializer::MemoryManager& memoryManager,
                              const SimulationSettings& settings) {
  SCOREP_USER_REGION("addClusters", SCOREP_USER_REGION_TYPE_FUNCTION);
  std::vector<std::unique_ptr<GhostCluster>> ghostClusters;

  // store the time stepping
  this->clusterLayout_ = clusterLayout;

  auto clusteringWriter = writer::ClusteringWriter(seissolInstance_.parameters().output.prefix);

  std::vector<std::size_t> drCellsPerCluster(clusterLayout.globalClusterCount);

  // setup DR schedulers
  for (const auto& layer : memoryManager.drStorage().leaves()) {
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

  const auto drOutputTimestep = drClusterOutput == std::numeric_limits<std::size_t>::max()
                                    ? std::numeric_limits<double>::infinity()
                                    : clusterLayout.timestepRate(drClusterOutput);

  std::vector<AbstractTimeCluster*> cellClusterBackmap(
      memoryManager.ltsStorage().getColorMap().size());

  const auto haloId = [&](const auto& id, HaloType halo) {
    auto cloned = id;
    cloned.halo = halo;
    return memoryManager.ltsStorage().getColorMap().colorId(cloned);
  };

  // iterate over local time clusters
  for (auto& layer : memoryManager.ltsStorage().leaves(Ghost)) {
    auto globalData = memoryManager.globalData();

    const auto clusterId = layer.getIdentifier().lts;

    // chop off at synchronization time
    const auto timeStepSize = clusterLayout.timestepRate(clusterId);
    const auto timeStepRate = clusterLayout.clusterRate(clusterId);

    // We print progress only if it is the cluster with the largest time step on each rank.
    // This does not mean that it is the largest cluster globally!
    const bool printProgress = (clusterId == clusterLayout.globalClusterCount - 1) &&
                               (layer.getIdentifier().halo == HaloType::Interior);
    const auto profilingId = layer.id();

    auto& cluster = cellClusters_.emplace_back(
        std::make_unique<CellCluster>(clusterId,
                                      clusterId,
                                      profilingId,
                                      settings,
                                      layer.getIdentifier().halo,
                                      timeStepSize,
                                      timeStepRate,
                                      printProgress,
                                      globalData,
                                      &layer,
                                      seissolInstance_,
                                      &loopStatistics_,
                                      &actorStateStatisticsManager_.addCluster(profilingId)));

    const auto clusterSize = layer.size();
    const auto dynRupSize = memoryManager.drStorage().layer(layer.id()).size();
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

  for (auto& layer : memoryManager.ltsStorage().leaves(Ghost)) {
    for (auto& otherLayer : memoryManager.ltsStorage().leaves(Ghost)) {
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

  // The dynamic rupture faces are laid out along the same colors as the cells, and both sides of a
  // fault face lie in the same time cluster. The faces of an interior layer lie between interior or
  // copy cells; the faces of a copy layer between copy and ghost cells.
  std::vector<DynamicRuptureCluster*> faceClusterBackmap(
      memoryManager.drStorage().getColorMap().size(), nullptr);
  for (auto& layer : memoryManager.drStorage().leaves(Ghost)) {
    if (layer.size() == 0) {
      continue;
    }

    const auto clusterId = layer.getIdentifier().lts;
    const auto profilingId = layer.id();
    auto* cellCluster = cellClusterBackmap[layer.id()];

    auto& cluster = faceClusters_.emplace_back(std::make_unique<DynamicRuptureCluster>(
        clusterLayout.timestepRate(clusterId),
        clusterLayout.clusterRate(clusterId),
        cellCluster->getExecutor(),
        profilingId,
        drOutputTimestep,
        memoryManager.globalData(),
        &layer,
        memoryManager.frictionLaw(),
        memoryManager.frictionLawDevice(),
        memoryManager.faultOutputManager(),
        seissolInstance_,
        &loopStatistics_,
        &actorStateStatisticsManager_.addCluster(profilingId)));

    cluster->setPriority(cellCluster->getPriority());

    cluster->connect(*cellClusterBackmap[haloId(layer.getIdentifier(), HaloType::Copy)]);
    if (layer.getIdentifier().halo == HaloType::Interior) {
      cluster->connect(*cellClusterBackmap[haloId(layer.getIdentifier(), HaloType::Interior)]);
    }

    faceClusterBackmap[layer.id()] = cluster.get();
  }

  // Create ghost time clusters for MPI
  followPlan_ = useTimeSteppingPlan(seissolInstance_.env());
  if (followPlan_) {
    logInfo() << "The clusters take their steps along the time stepping plan.";
  }

  haloTransports_ =
      std::make_unique<HaloTransportFactory>(Mpi::mpi.getPreferredDataTransferMode(),
                                             usePersistentMpi(seissolInstance_.env()),
                                             useCclPerDirection(seissolInstance_.env()),
                                             clusterLayout.globalClusterCount);
  for (auto& layer : memoryManager.ltsStorage().leaves(Ghost | Interior)) {

    const auto displayName = "copy-" + std::to_string(layer.getIdentifier().lts);

    for (const auto [i, halo] : common::enumerate(haloStructure.at(layer.id()))) {

      const bool hasNeighborRegions = !halo.copy.empty() || !halo.ghost.empty();
      const auto other = memoryManager.ltsStorage().getColorMap().argument(i);

      if (hasNeighborRegions) {

        assert(other.halo == HaloType::Ghost);
        assert(other.lts + 1 >= layer.getIdentifier().lts);
        assert(other.lts <
               std::min(layer.getIdentifier().lts + 2, clusterLayout.globalClusterCount));

        const auto otherTimeStepSize = clusterLayout.timestepRate(other.lts);
        const auto otherTimeStepRate = clusterLayout.clusterRate(other.lts);

        const auto otherDisplayName = "ghost-" + std::to_string(other.lts);

        const auto& regions = haloStructure.at(layer.id()).at(i);
        ghostClusters.push_back(std::make_unique<GhostCluster>(
            otherTimeStepSize,
            otherTimeStepRate,
            displayName,
            otherDisplayName,
            regions,
            haloTransports_->create(regions, layer.getIdentifier().lts, other.lts)));

        // Connect with previous copy layer.
        ghostClusters.back()->connect(*cellClusterBackmap[layer.id()]);

        // the dynamic rupture faces of the copy layer read the ghost cells of their own cluster
        auto* faceCluster = faceClusterBackmap[layer.id()];
        if (faceCluster != nullptr && other.lts == layer.getIdentifier().lts) {
          faceCluster->observe(*ghostClusters.back());
        }
      }
    }
  }

  clusteringWriter.write();

  reportClusterDependencies(memoryManager, haloStructure);

  // Sort clusters by time step size in increasing order
  auto rateSorter = [](const auto& a, const auto& b) {
    return a->getTimeStepRate() < b->getTimeStepRate();
  };
  std::sort(cellClusters_.begin(), cellClusters_.end(), rateSorter);

  for (const auto& cluster : cellClusters_) {
    clusters_.emplace_back(cluster.get());
  }
  for (const auto& cluster : faceClusters_) {
    clusters_.emplace_back(cluster.get());
  }
  std::stable_sort(clusters_.begin(), clusters_.end(), rateSorter);

  for (auto* cluster : clusters_) {
    if (cluster->getPriority() == ActorPriority::High) {
      highPrioClusters_.emplace_back(cluster);
    } else {
      lowPrioClusters_.emplace_back(cluster);
    }
  }

  std::sort(ghostClusters.begin(), ghostClusters.end(), rateSorter);

  if (seissol::useCommThread(Mpi::mpi, seissolInstance_.env())) {
    communicationManager_ = std::make_unique<ThreadedCommunicationManager>(
        std::move(ghostClusters), &seissolInstance_.pinning());
  } else {
    communicationManager_ = std::make_unique<SerialCommunicationManager>(std::move(ghostClusters));
  }

  auto* ghostClusterPointer = communicationManager_->getGhostClusters();

  concurrent_ = isDeviceOn() && useConcurrentClusters(seissolInstance_.env());
  if (concurrent_) {
    logInfo() << "The clusters run concurrently on the device.";
    for (auto* cluster : clusters_) {
      cluster->setConcurrent(true);
    }
    for (auto& cluster : *ghostClusterPointer) {
      cluster->setConcurrent(true);
    }
  }

  std::vector<AbstractTimeCluster*> allClusters(clusters_.size() + ghostClusterPointer->size());
  for (std::size_t i = 0; i < clusters_.size(); ++i) {
    allClusters[i] = clusters_[i];
  }
  for (std::size_t i = 0; i < ghostClusterPointer->size(); ++i) {
    allClusters[i + clusters_.size()] = ghostClusterPointer->at(i).get();
  }

  for (auto& timeMirrorManager : seissolInstance_.getTimeMirrorManagers()) {
    timeMirrorManager.setClusterVector(allClusters);
  }
}

void TimeManager::setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
  this->faultOutputManager_ = faultOutputManager;
  for (auto& cluster : faceClusters_) {
    cluster->setFaultOutputManager(faultOutputManager);
  }
}

seissol::dr::output::OutputManager* TimeManager::faultOutputManager() {
  assert(faultOutputManager_ != nullptr);
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

  if (followPlan_) {
    followPlan();
  }

  // Also completes the synchronization after following the plan.
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
  if (concurrent_) {
    // the clusters have only enqueued their work
    device.api->syncDevice();
  }
  device.api->popLastProfilingMark();
#endif
  for (auto& cluster : clusters_) {
    cluster->finishPhase();
  }
}

void TimeManager::followPlan() {
  std::vector<PlannedCluster> planned;
  planned.reserve(clusters_.size());
  for (auto* cluster : clusters_) {
    planned.push_back({cluster->getTimeStepRate(),
                       cluster->getStepsUntilSync(),
                       cluster->dataReadiness(),
                       cluster->getPriority()});
  }

  long largestRate = 1;
  long lastTick = 0;
  for (const auto& cluster : planned) {
    largestRate = std::max(largestRate, cluster.timeStepRate);
    lastTick = std::max(lastTick, cluster.stepsUntilSync);
  }
  // per super-timestep: whether it takes output samples, and whether it has any irregular step
  std::map<long, std::pair<bool, bool>> superStepWork;

  for (const auto& step : planTimeSteps(planned)) {
    auto* cluster = clusters_[step.cluster];
    // along the plan, a cluster can only have to wait for the halo exchange
    auto action = cluster->getNextLegalAction();
    while (action == ActorAction::Nothing) {
      communicationManager_->progression();
      action = cluster->getNextLegalAction();
    }
    if (action != step.action) {
      logError() << "The cluster" << cluster->identifier() << "is not ready for step" << step.step
                 << "of the time stepping plan.";
    }
    cluster->act();

    // the super-timestep the action falls into, by the first tick it covers
    const auto rate = planned[step.cluster].timeStepRate;
    const auto superStep = step.step * rate / largestRate;
    auto& [outputs, irregular] = superStepWork[superStep];
    outputs = outputs || cluster->lastStepWork().outputs;
    irregular = irregular || cluster->lastStepWork().irregular();
  }

  for (const auto& [superStep, work] : superStepWork) {
    ++superSteps_;
    if ((superStep + 1) * largestRate > lastTick) {
      ++shortenedSuperSteps_;
    } else {
      outputFreeSuperSteps_ += work.first ? 0 : 1;
      regularSuperSteps_ += work.second ? 0 : 1;
    }
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
  for (auto& cluster : cellClusters_) {
    cluster->setPointSources(std::move(sourceClusters[cluster->layerId()]));
  }
}

void TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter) {
  for (auto& cluster : cellClusters_) {
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

  // every message sent has to be received somewhere
  std::array<std::size_t, 2> messages{0, 0};
  for (auto& cluster : *communicationManager_->getGhostClusters()) {
    messages[0] += cluster->sentMessages();
    messages[1] += cluster->receivedMessages();
    cluster->finalize();
  }
  MPI_Allreduce(MPI_IN_PLACE,
                messages.data(),
                messages.size(),
                Mpi::castToMpiType<std::size_t>(),
                MPI_SUM,
                Mpi::mpi.comm());
  logInfo() << "Halo exchange:" << messages[0] << "messages sent," << messages[1]
            << "received (summed over all ranks)";

  if (followPlan_) {
    std::array<std::size_t, 4> superSteps{
        superSteps_, shortenedSuperSteps_, outputFreeSuperSteps_, regularSuperSteps_};
    MPI_Allreduce(MPI_IN_PLACE,
                  superSteps.data(),
                  superSteps.size(),
                  Mpi::castToMpiType<std::size_t>(),
                  MPI_SUM,
                  Mpi::mpi.comm());
    logInfo() << "Super-timesteps:" << superSteps[0]
              << "; ending at a synchronization point:" << superSteps[1]
              << "; full ones without output samples:" << superSteps[2]
              << "; full ones without output samples or host work:" << superSteps[3]
              << "(summed over all ranks)";
  }
  if (messages[0] != messages[1]) {
    logWarning() << "The halo exchange sent and received a different number of messages.";
  }

  communicationManager_.reset(nullptr);
  haloTransports_.reset(nullptr);
}

void TimeManager::synchronizeTo(seissol::initializer::AllocationPlace place) {
#ifdef ACL_DEVICE
  bool sameExecutor = true;
  for (auto& cluster : clusters_) {
    sameExecutor &= clusters_.front()->getExecutor() == cluster->getExecutor();
  }
  if (sameExecutor) {
    seissolInstance_.memoryManager().synchronizeTo(place);
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

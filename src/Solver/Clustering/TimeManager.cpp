// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 *http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 **/

#include "Parallel/MPI.h"

#include "Parallel/Helper.h"
#include "ResultWriter/ClusteringWriter.h"
#include "SeisSol.h"
#include "Solver/Clustering/Communication/CommunicationManager.h"
#include "TimeManager.h"
#include <AbstractAPI.h>
#include <Common/Executor.h>
#include <DynamicRupture/Output/OutputManager.h>
#include <Initializer/MemoryManager.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Tree/Layer.h>
#include <Initializer/Typedefs.h>
#include <Kernels/PointSourceCluster.h>
#include <Parallel/Host/SyncExecutor.h>
#include <ResultWriter/ReceiverWriter.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <Solver/Clustering/Communication/CommunicationFactory.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <Solver/Clustering/Computation/CopyCluster.h>
#include <Solver/Clustering/Computation/DynamicRuptureCluster.h>
#include <Solver/Clustering/Computation/GhostCluster.h>
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

namespace seissol::solver::clustering {

TimeManager::TimeManager(seissol::SeisSol& seissolInstance)
    : logUpdates(std::numeric_limits<unsigned int>::max()), seissolInstance(seissolInstance),
      actorStateStatisticsManager(loopStatistics) {
  loopStatistics.addRegion("computeLocalIntegration");
  loopStatistics.addRegion("computeNeighboringIntegration");
  loopStatistics.addRegion("computeDynamicRupture");
  loopStatistics.addRegion("computePointSources");

  loopStatistics.enableSampleOutput(
      seissolInstance.getSeisSolParameters().output.loopStatisticsNetcdfOutput);

  cpuExecutor = std::make_shared<parallel::host::SyncExecutor>();
}

void TimeManager::addClusters(initializer::ClusterLayout& layout,
                              const communication::HaloCommunication& halo,
                              initializer::MemoryManager& memoryManager,
                              bool usePlasticity) {
  const auto globalData = memoryManager.getGlobalData();
  int profilingId = 0;

  this->layout = layout;

  auto clusteringWriter = writer::ClusteringWriter(memoryManager.getOutputPrefix());

  communication::CommunicationClusterFactory communicationFactory(
      communication::CommunicationMode::DirectMPI, layout.globalClusterCount);
  communicationFactory.prepare();

  const auto sendClusters = communicationFactory.getAllSends(halo, cpuExecutor, PriorityHighest);
  const auto recvClusters = communicationFactory.getAllRecvs(halo, cpuExecutor, PriorityHighest);

  std::vector<std::unordered_map<LayerType, std::shared_ptr<AbstractTimeCluster>>>
      cellClusterBackmap(layout.localClusterIds.size());
  std::vector<std::unordered_map<LayerType, std::shared_ptr<AbstractTimeCluster>>>
      faceClusterBackmap(layout.localClusterIds.size());

  for (auto& layer : memoryManager.getLtsTree()->leaves()) {
    if (layer.getNumberOfCells() == 0) {
      continue;
    }

    const auto globalClusterId = layer.getClusterId();
    const auto localClusterId = layout.localClusterFromGlobal(globalClusterId).value();
    const auto timeStepRate = layout.clusterRate(localClusterId);
    const auto timeStepSize = layout.timestepRate(localClusterId);
    if (layer.getLayerType() == Interior) {
      const bool printProgress = localClusterId == *std::max_element(layout.localClusterIds.begin(),
                                                                     layout.localClusterIds.end());
      // add interior cluster
      clusters.push_back(std::make_shared<computation::TimeCluster>(
          localClusterId,
          globalClusterId,
          profilingId,
          usePlasticity,
          timeStepSize,
          timeStepRate,
          printProgress,
          globalData,
          &layer,
          memoryManager.getLts(),
          seissolInstance,
          &loopStatistics,
          &actorStateStatisticsManager.addCluster(profilingId),
          cpuExecutor,
          PriorityLowest));
      ++profilingId;
      lowPrioClusters.push_back(clusters.back());
    }
    if (layer.getLayerType() == Copy) {
      clusters.push_back(std::make_shared<computation::CopyCluster>(
          localClusterId,
          globalClusterId,
          profilingId,
          usePlasticity,
          timeStepSize,
          timeStepRate,
          globalData,
          &layer,
          memoryManager.getLts(),
          seissolInstance,
          &loopStatistics,
          &actorStateStatisticsManager.addCluster(profilingId),
          sendClusters.at(globalClusterId),
          cpuExecutor,
          PriorityNormal));
      ++profilingId;
      highPrioClusters.push_back(clusters.back());
    }
    if (layer.getLayerType() == Ghost) {
      clusters.push_back(
          std::make_shared<computation::GhostCluster>(timeStepSize,
                                                      timeStepRate,
                                                      recvClusters.at(globalClusterId),
                                                      cpuExecutor,
                                                      PriorityNormal));
      highPrioClusters.push_back(clusters.back());
    }

    cellClusterBackmap[localClusterId][layer.getLayerType()] = clusters.back();
    if (layer.getLayerType() != Ghost) {
      cellClusters.push_back(std::dynamic_pointer_cast<computation::TimeCluster>(clusters.back()));

      clusteringWriter.addCluster(
          profilingId, localClusterId, layer.getLayerType(), layer.getNumberOfCells(), 0);
    }
  }
  for (auto& layer : memoryManager.getDynamicRuptureTree()->leaves(Ghost)) {
    if (layer.getNumberOfCells() == 0) {
      continue;
    }

    const auto globalClusterId = layer.getClusterId();
    const auto localClusterId = layout.localClusterFromGlobal(globalClusterId).value();
    const auto timeStepRate = layout.clusterRate(localClusterId);
    const auto timeStepSize = layout.timestepRate(localClusterId);
    if (layer.getLayerType() == Interior || layer.getLayerType() == Copy) {
      clusters.push_back(std::make_shared<computation::DynamicRuptureCluster>(
          timeStepSize,
          timeStepRate,
          profilingId,
          layer.getLayerType(),
          &layer,
          memoryManager.getDynamicRupture(),
          globalData,
          memoryManager.getFrictionLaw(),
          memoryManager.getFrictionLawDevice(),
          memoryManager.getFaultOutputManager(),
          seissolInstance,
          &loopStatistics,
          &actorStateStatisticsManager.addCluster(profilingId),
          cpuExecutor,
          layer.getLayerType() == Interior ? PriorityLowest : PriorityNormal));
      ++profilingId;
      if (layer.getLayerType() == Interior) {
        lowPrioClusters.push_back(clusters.back());
      } else {
        highPrioClusters.push_back(clusters.back());
      }
    }

    faceClusterBackmap[localClusterId][layer.getLayerType()] = clusters.back();
    faceClusters.push_back(
        std::dynamic_pointer_cast<computation::DynamicRuptureCluster>(clusters.back()));

    clusteringWriter.addCluster(
        profilingId, localClusterId, layer.getLayerType(), layer.getNumberOfCells(), 0);
  }

  const auto connectIfBothExist = [](auto& a, auto& b) {
    if (a && b) {
      a->connect(*b);
    }
  };

  for (std::size_t i = 0; i < layout.localClusterIds.size(); ++i) {
    connectIfBothExist(cellClusterBackmap[i][Copy], cellClusterBackmap[i][Interior]);
    connectIfBothExist(cellClusterBackmap[i][Copy], cellClusterBackmap[i][Ghost]);

    connectIfBothExist(faceClusterBackmap[i][Interior], cellClusterBackmap[i][Interior]);
    connectIfBothExist(faceClusterBackmap[i][Interior], cellClusterBackmap[i][Copy]);

    connectIfBothExist(faceClusterBackmap[i][Copy], cellClusterBackmap[i][Copy]);
    connectIfBothExist(faceClusterBackmap[i][Copy], cellClusterBackmap[i][Ghost]);
  }

  for (std::size_t i = 1; i < layout.localClusterIds.size(); ++i) {
    if (layout.localClusterIds[i - 1] + 1 == layout.localClusterIds[i]) {
      connectIfBothExist(cellClusterBackmap[i][Interior], cellClusterBackmap[i - 1][Interior]);
      connectIfBothExist(cellClusterBackmap[i][Interior], cellClusterBackmap[i - 1][Copy]);
      connectIfBothExist(cellClusterBackmap[i][Copy], cellClusterBackmap[i - 1][Interior]);
      connectIfBothExist(cellClusterBackmap[i][Copy], cellClusterBackmap[i - 1][Copy]);
    }
  }

  std::vector<std::shared_ptr<clustering::communication::NeighborCluster>> communication;

  communication.insert(communication.end(), sendClusters.begin(), sendClusters.end());
  communication.insert(communication.end(), recvClusters.begin(), recvClusters.end());

  if (seissol::useCommThread(MPI::mpi)) {
    communicationManager = std::make_shared<communication::ThreadedCommunicationManager>(
        communication, &seissolInstance.getPinning());
  } else {
    communicationManager =
        std::make_shared<communication::SerialCommunicationManager>(communication);
  }

  auto& timeMirrorManagers = seissolInstance.getTimeMirrorManagers();
  auto& [increaseManager, decreaseManager] = timeMirrorManagers;

  increaseManager.setTimeClusterVector(clusters);
  decreaseManager.setTimeClusterVector(clusters);
}

void TimeManager::setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
  this->faultOutputManager = faultOutputManager;
  for (auto& cluster : faceClusters) {
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

  communicationManager->reset();

  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  device.api->putProfilingMark("advanceInTime", device::ProfilingColors::Blue);
#endif

  const auto epochExecute = [&]() {
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
      for (int i = 0; i < 5; ++i) {
        std::for_each(highPrioClusters.begin(), highPrioClusters.end(), [&](auto& cluster) {
          communicationManager->progression();
          cluster->act();
        });
      }
      std::for_each(lowPrioClusters.begin(), lowPrioClusters.end(), [&](auto& cluster) {
        communicationManager->progression();
        cluster->act();
      });
      finished =
          std::all_of(clusters.begin(), clusters.end(), [](auto& c) { return c->synchronized(); });
      finished &= communicationManager->checkIfFinished();
    }
  };

  cpuExecutor->start([&](const auto&) { epochExecute(); }, &seissolInstance.getPinning());

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

double TimeManager::getTimeTolerance() { return 1E-5 * layout.value().minimumTimestep; }

void TimeManager::setPointSourcesForClusters(
    std::unordered_map<LayerType, std::vector<seissol::kernels::PointSourceClusterPair>>
        sourceClusters) {
  for (auto& cluster : cellClusters) {
    auto layerClusters = sourceClusters.find(cluster->getLayerType());
    if (layerClusters != sourceClusters.end() &&
        cluster->getClusterId() < layerClusters->second.size()) {
      cluster->setPointSources(std::move(layerClusters->second.at(cluster->getClusterId())));
    }
  }
}

void TimeManager::setReceiverClusters(writer::ReceiverWriter& receiverWriter) {
  for (auto& cluster : cellClusters) {
    cluster->setReceiverCluster(
        receiverWriter.receiverCluster(cluster->getClusterId(), cluster->getLayerType()));
  }
}

void TimeManager::setInitialTimes(double time) {
  assert(time >= 0);

  for (auto& cluster : clusters) {
    cluster->setTime(time);
  }

  for (auto& cluster : cellClusters) {
    cluster->setReceiverTime(time);
  }
}

void TimeManager::setTv(double tv) {
  for (auto& cluster : cellClusters) {
    cluster->setTv(tv);
  }
}

void TimeManager::freeDynamicResources() {
  for (auto& cluster : clusters) {
    cluster->finalize();
  }
  communicationManager.reset();
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
    device::DeviceInstance::getInstance().api->syncDefaultStreamWithHost();
  }
#endif
}

} // namespace seissol::solver::clustering

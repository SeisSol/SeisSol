// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "TimeManager.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Executor.h"
#include "Common/Iterator.h"
#include "Common/Json.h"
#include "Common/Typedefs.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "DynamicRupture/Output/PickpointSampling.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Kernels/Common.h"
#include "Kernels/PointSourceCluster.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/HostData.h"
#include "Memory/HostSnapshot.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/DeviceTimeline.h"
#include "Monitoring/DeviceTimer.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/RunSummary.h"
#include "Monitoring/Trace.h"
#include "Monitoring/Unit.h"
#include "Monitoring/Watchdog.h"
#include "Numerical/Statistics.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/OpenMP.h"
#include "Parallel/Pin.h"
#include "Parallel/Runtime/HostActor.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Parallel/Runtime/Stream.h"
#include "ResultWriter/ClusteringWriter.h"
#include "ResultWriter/ReceiverWriter.h"
#include "SeisSol.h"
#include "Solver/Settings.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Compute/CellCluster.h"
#include "Solver/TimeStepping/Compute/DynamicRuptureCluster.h"
#include "Solver/TimeStepping/Halo/CommunicationManager.h"
#include "Solver/TimeStepping/Halo/GhostCluster.h"
#include "Solver/TimeStepping/Halo/HaloCheck.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/TaskHaloTransport.h"
#include "Solver/TimeStepping/HaloTransportFactory.h"
#include "Solver/TimeStepping/Plan/SuperStepRecorder.h"
#include "Solver/TimeStepping/Plan/TimeSteppingPlan.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <pthread.h>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/AbstractAPI.h>
#include <Device/device.h>
#endif

namespace seissol::solver {

TimeManager::TimeManager(seissol::SeisSol& seissolInstance)
    : seissolInstance_(seissolInstance), actorStateStatisticsManager_(loopStatistics_) {
  loopStatistics_.addRegion("computeLocalIntegration");
  loopStatistics_.addRegion("computeNeighboringIntegration");
  loopStatistics_.addRegion("computeDynamicRupture");
  loopStatistics_.addRegion("computePointSources");
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
                               const HaloCommunication& haloStructure) {
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

/**
 * By rank, whether the process is on the own node. For tests, SEISSOL_EXCHANGE_RANKS_PER_NODE=N
 * makes nodes of N consecutive ranks each instead.
 */
std::vector<bool> ranksOnNode(utils::Env& env) {
  const auto size = static_cast<std::size_t>(Mpi::mpi.size());
  const auto rank = static_cast<std::size_t>(Mpi::mpi.rank());
  std::vector<bool> onNode(size, false);
  const auto perNode = env.get<std::size_t>("EXCHANGE_RANKS_PER_NODE", 0);
  if (perNode > 0) {
    for (std::size_t other = 0; other < size; ++other) {
      onNode[other] = other / perNode == rank / perNode;
    }
    return onNode;
  }
  MPI_Group world{};
  MPI_Group shared{};
  MPI_Comm_group(Mpi::mpi.comm(), &world);
  MPI_Comm_group(Mpi::mpi.sharedMemComm(), &shared);
  int sharedSize = 0;
  MPI_Group_size(shared, &sharedSize);
  std::vector<int> sharedRanks(static_cast<std::size_t>(sharedSize));
  std::iota(sharedRanks.begin(), sharedRanks.end(), 0);
  std::vector<int> worldRanks(sharedRanks.size());
  MPI_Group_translate_ranks(shared, sharedSize, sharedRanks.data(), world, worldRanks.data());
  MPI_Group_free(&shared);
  MPI_Group_free(&world);
  for (const auto other : worldRanks) {
    onNode.at(static_cast<std::size_t>(other)) = true;
  }
  return onNode;
}

} // namespace

void TimeManager::addClusters(const initializer::ClusterLayout& clusterLayout,
                              const HaloCommunication& haloStructure,
                              initializer::MemoryManager& memoryManager,
                              const SimulationSettings& settings) {
  SCOREP_USER_REGION("addClusters", SCOREP_USER_REGION_TYPE_FUNCTION);
  std::vector<std::unique_ptr<GhostCluster>> ghostClusters;

  setUpHostScheduling();
  DeviceTimer::available() = DeviceTimer::canTimeEvents();
  DeviceTimer::enabled() = useDeviceLoopStatistics(seissolInstance_.env());

  // store the time stepping
  this->clusterLayout_ = clusterLayout;

  // written in initIO, with the rest of the output
  auto& clusteringWriter = clusteringWriter_.emplace(seissolInstance_.parameters().output.prefix,
                                                     seissolInstance_.runId());

  std::vector<std::size_t> drCellsPerCluster(clusterLayout.globalClusterCount);

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

  // the on-fault receivers sample every pickdt, or else every printtimeinterval time steps of the
  // finest cluster with fault faces, counted from time zero
  const auto& pickpointParameters = seissolInstance_.parameters().output.pickpointParameters;
  const double pickpointInterval = pickpointParameters.samplingInterval > 0
                                       ? pickpointParameters.samplingInterval
                                       : pickpointParameters.printTimeInterval * drOutputTimestep;
  const dr::output::PickpointSampling pickpointSampling(
      pickpointInterval, seissolInstance_.parameters().timeStepping.endTime);

  std::vector<AbstractTimeCluster*> cellClusterBackmap(
      memoryManager.ltsStorage().getColorMap().size());

  const auto haloId = [&](const auto& id, HaloType halo) {
    auto cloned = id;
    cloned.halo = halo;
    return memoryManager.ltsStorage().getColorMap().colorId(cloned);
  };

  // iterate over local time clusters
  for (auto& layer : memoryManager.ltsStorage().leaves(Ghost)) {
    const auto clusterId = layer.getIdentifier().lts;

    // chop off at synchronization time
    const auto timeStepSize = clusterLayout.timestepRate(clusterId);
    const auto timeStepRate = clusterLayout.clusterRate(clusterId);

    // We print progress only if it is the cluster with the largest time step on each rank.
    // This does not mean that it is the largest cluster globally!
    const bool printProgress = (clusterId == clusterLayout.globalClusterCount - 1) &&
                               (layer.getIdentifier().halo == HaloType::Interior);
    const auto profilingId = layer.id();

    // the cluster computes in the configuration of its layer
    auto& cluster = cellClusters_.emplace_back(dispatchConfig(
        layer.getIdentifier().config, [&](auto cfg) -> std::unique_ptr<CellClusterInterface> {
          using Cfg = decltype(cfg);
          return std::make_unique<CellCluster<Cfg>>(
              clusterId,
              clusterId,
              profilingId,
              settings,
              layer.getIdentifier().halo,
              timeStepSize,
              timeStepRate,
              printProgress,
              memoryManager.globalData<Cfg>(),
              &layer,
              seissolInstance_,
              &loopStatistics_,
              &deviceTimeline_,
              &actorStateStatisticsManager_.addCluster(profilingId));
        }));

    const auto clusterSize = layer.size();
    const auto dynRupSize = memoryManager.drStorage().layer(layer.id()).size();
    // a copy cell is stored once per neighbor it sends to, and every copy is computed
    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    std::size_t duplicateCells = 0;
    for (std::size_t cell = 0; cell < clusterSize; ++cell) {
      if (secondaryInformation[cell].duplicate != 0) {
        ++duplicateCells;
      }
    }
    // the dynamic rupture faces of the layer get counted with their face cluster
    clusterSizes_.push_back(
        ClusterSize{cluster.get(), clusterId, clusterSize - duplicateCells, duplicateCells, 0});
    // Add writer to output
    clusteringWriter.addCluster(profilingId,
                                clusterId,
                                layer.getIdentifier().halo,
                                clusterSize,
                                dynRupSize,
                                timeStepSize,
                                static_cast<long>(timeStepRate));

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
  std::vector<DynamicRuptureClusterInterface*> faceClusterBackmap(
      memoryManager.drStorage().getColorMap().size(), nullptr);
  for (auto& layer : memoryManager.drStorage().leaves(Ghost)) {
    if (layer.size() == 0) {
      continue;
    }

    const auto clusterId = layer.getIdentifier().lts;
    const auto profilingId = layer.id();
    auto* cellCluster = cellClusterBackmap[layer.id()];

    auto& cluster = faceClusters_.emplace_back(
        dispatchConfig(layer.getIdentifier().config,
                       [&](auto cfg) -> std::unique_ptr<DynamicRuptureClusterInterface> {
                         using Cfg = decltype(cfg);
                         return std::make_unique<DynamicRuptureCluster<Cfg>>(
                             clusterLayout.timestepRate(clusterId),
                             clusterLayout.clusterRate(clusterId),
                             cellCluster->getExecutor(),
                             profilingId,
                             pickpointSampling,
                             memoryManager.globalData<Cfg>(),
                             &layer,
                             memoryManager.frictionLaw(),
                             memoryManager.frictionLawDevice(),
                             memoryManager.faultOutputManager(),
                             seissolInstance_,
                             &loopStatistics_,
                             &deviceTimeline_,
                             &actorStateStatisticsManager_.addCluster(profilingId));
                       }));

    cluster->setPriority(cellCluster->getPriority());
    clusterSizes_.push_back(ClusterSize{cluster.get(), clusterId, 0, 0, layer.size()});

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

  // the transports go between the time clusters of each configuration of the run; the schedulers
  // of the device streams tell their directions apart by both, the same way on all processes
  const auto runConfigs = seissolInstance_.parameters().model.configs();
  const auto haloCluster = [&](const initializer::LayerIdentifier& identifier) {
    const auto position = static_cast<std::size_t>(
        std::find(runConfigs.begin(), runConfigs.end(), identifier.config) - runConfigs.begin());
    assert(position < runConfigs.size());
    return identifier.lts * runConfigs.size() + position;
  };
  haloTransports_ =
      std::make_unique<HaloTransportFactory>(Mpi::mpi.getPreferredDataTransferMode(),
                                             usePersistentMpi(seissolInstance_.env()),
                                             useExchangePerDirection(seissolInstance_.env()),
                                             clusterLayout.globalClusterCount * runConfigs.size(),
                                             runConfigs.size(),
                                             exchangeLanes(seissolInstance_.env()));
  if (Mpi::mpi.getPreferredDataTransferMode() == Mpi::DataTransferMode::ActorMpi) {
    // MPI goes on while this thread waits for the device; on the CPUs the OpenMP threads leave
    // free, if there are any
    const auto cpus = seissolInstance_.pinning().getFreeCPUsMask();
    std::function<void()> pin;
    if (!parallel::Pinning::freeCPUsMaskEmpty(cpus)) {
      pin = [cpus]() { cpus.setaffinity_np(pthread_self()); };
    }
    parallel::runtime::HostActor::instance().startProgressThread(pin);
    logInfo() << "A thread of its own" << (pin ? "on the free CPUs" : "")
              << "sends and receives the halo data with MPI through host memory; the device"
              << "waits for it.";
  }
  if (useExchangeSplitByNode(seissolInstance_.env()) && haloTransports_->scheduler() != nullptr) {
    haloTransports_->splitByNode(ranksOnNode(seissolInstance_.env()));
    logInfo() << "The halo exchange on device streams stays within the node; with the processes"
              << "on other nodes, it goes over GPU-aware MPI from the host.";
  }
  // what the test of the halo exchange before the time stepping goes through
  std::vector<CheckedTransport> checkedTransports;
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
        auto transport = haloTransports_->create(
            regions, haloCluster(layer.getIdentifier()), haloCluster(other));
        auto description = displayName;
        description.append(" and ").append(otherDisplayName);
        checkedTransports.push_back(
            {&transport->dataTransport(), &regions, std::move(description)});
        ghostClusters.push_back(std::make_unique<GhostCluster>(otherTimeStepSize,
                                                               otherTimeStepRate,
                                                               displayName,
                                                               otherDisplayName,
                                                               regions,
                                                               std::move(transport)));

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

  // the library of the halo exchange knows all transports now; collective over all processes
  if (haloTransports_->scheduler() != nullptr) {
    haloTransports_->scheduler()->prepare();
  }
  if (useHaloCheck(seissolInstance_.env())) {
    // before the exchanges get timed, which would count the test
    const auto watchdog = monitoring::Watchdog::settingsFromEnv(seissolInstance_.env());
    testHaloExchange(checkedTransports,
                     isDeviceOn() ? HaloMemory::Device : HaloMemory::Host,
                     Mpi::dataTransferModeName(Mpi::mpi.getPreferredDataTransferMode()),
                     watchdog.enabled ? watchdog.timeout : std::chrono::seconds(0),
                     watchdog.abort);
  }
  if (haloTransports_->scheduler() != nullptr && DeviceTimeline::active()) {
    exchangeTimer_ =
        std::make_unique<DeviceTimer>(nullptr, &deviceTimeline_, DeviceTimeline::Kind::Exchanges);
    exchangeTimer_->setTrace([] { return std::string("halo exchange"); }, "exchange");
    haloTransports_->scheduler()->setTimer(exchangeTimer_.get(), &exchangeDeviceTime_);
  }

  commThread_ = seissol::useCommThread(Mpi::mpi, seissolInstance_.env());
  if (haloTransports_->taskEngine() != nullptr) {
    // the main thread creates the tasks of the exchanges; where the OpenMP runtime may wait for one
    // of them right there, it cannot complete them as well
    const bool pollingThread =
        commThread_ || !parallel::runtime::HostRuntime::detachedTasksNeverBlock();
    if (!commThread_ && pollingThread && !ghostClusters.empty()) {
      logInfo() << "A thread of its own polls the halo exchange in the task graph, as this OpenMP"
                << "runtime may wait for the tasks of the exchanges where they get created.";
    }
    communicationManager_ =
        std::make_unique<TaskCommunicationManager>(std::move(ghostClusters),
                                                   *haloTransports_->taskEngine(),
                                                   &seissolInstance_.pinning(),
                                                   pollingThread);
  } else if (commThread_) {
    communicationManager_ = std::make_unique<ThreadedCommunicationManager>(
        std::move(ghostClusters), &seissolInstance_.pinning());
  } else {
    communicationManager_ = std::make_unique<SerialCommunicationManager>(std::move(ghostClusters));
  }

  auto* ghostClusterPointer = communicationManager_->getGhostClusters();

  if (isDeviceOn() && useHostActor(seissolInstance_.env())) {
    logInfo() << "The work on the host for the device (the samples of the receivers and fault"
              << "receivers) runs once the device is done with what it reads"
              << "(SEISSOL_HOST_ACTOR=1).";
  }

  concurrent_ = isDeviceOn() && useConcurrentClusters(seissolInstance_.env());
  if (concurrent_) {
    logInfo() << "The clusters run concurrently on the device.";
    parallel::runtime::hostFunctionsOnStreams() = true;
    for (auto* cluster : clusters_) {
      cluster->setConcurrent(true);
    }
    for (auto& cluster : *ghostClusterPointer) {
      cluster->setConcurrent(true);
    }
    auto* scheduler = haloTransports_->scheduler();
    if (scheduler != nullptr && scheduler->launchOrder() == LaunchOrder::Global) {
      scheduler->setStreamOrdered(true);
      logInfo() << "The halo exchange is ordered on the device.";
      if (useExchangeTransportEvents(seissolInstance_.env())) {
        scheduler->setTransportEvents(true);
        logInfo() << "The work of a cluster waits only for the halo exchange of its own transports"
                  << "(SEISSOL_EXCHANGE_TRANSPORT_EVENTS=1; with one lane only).";
      }
    }
  }

  if (useSuperStepGraphs(seissolInstance_.env())) {
    if (!followPlan_ || !concurrent_) {
      logWarning() << "Recording super-timesteps needs a device, SEISSOL_TIMESTEPPING_PLAN=1 and"
                   << "SEISSOL_CONCURRENT_CLUSTERS=1.";
    } else if (!ghostClusterPointer->empty() &&
               (haloTransports_->scheduler() == nullptr ||
                !haloTransports_->scheduler()->streamOrdered() || commThread_)) {
      logWarning() << "Recording super-timesteps with halo exchange needs a transfer mode on device"
                   << "streams (ccl, stream-mpi or shmem) in its global order, and no"
                   << "communication thread.";
    } else if (!ghostClusterPointer->empty() &&
               (!haloTransports_->scheduler()->recordable() || haloTransports_->split())) {
      logWarning() << "The halo exchange of this transfer mode cannot be replayed from a recording"
                   << "of super-timesteps (stream-mpi, shmem with Intel SHMEM, and any split by"
                   << "node); use ccl, or shmem with NVSHMEM or ROCSHMEM.";
    } else if (!SuperStepRecorder::available()) {
      logWarning() << "Recording super-timesteps needs a device that builds graphs node by node"
                   << "(and SEISSOL_DEVICE_GRAPH_NODES not set to 0).";
    } else {
      replay_ = true;
      const bool branches = useSuperStepBranches(seissolInstance_.env());
      recorder_ = std::make_unique<SuperStepRecorder>(branches);
      if (!branches) {
        logInfo() << "The branches of the work of a cluster run one after another in the recordings"
                  << "(SEISSOL_SUPERSTEP_BRANCHES=0).";
      }
      if (DeviceTimeline::active()) {
        superStepTimer_ = std::make_unique<DeviceTimer>(
            nullptr, &deviceTimeline_, DeviceTimeline::Kind::SuperSteps);
        superStepTimer_->setTrace([] { return std::string("super-timesteps"); }, "super-timestep");
        recorder_->setTimer(superStepTimer_.get());
      }
      graphMode_ = superStepGraphMode(seissolInstance_.env());
      if (graphMode_ == SuperStepGraphMode::Rebuild) {
        logInfo() << "Every super-timestep gets recorded anew and launched once, none replayed"
                  << "(SEISSOL_SUPERSTEP_GRAPH_MODE=rebuild).";
      } else {
        logInfo() << "Regular super-timesteps get recorded and replayed.";
      }
      if (graphMode_ == SuperStepGraphMode::Always) {
        logInfo() << "All other super-timesteps get recorded and launched once"
                  << "(SEISSOL_SUPERSTEP_GRAPH_MODE=always).";
      }

      // the outputs must not stand in the way of a recording
      for (auto* cluster : clusters_) {
        cluster->setRunTimeOutputs(true);
      }
      logInfo() << "The outputs decide about their samples when the work of a step runs.";
      if (std::any_of(clusters_.begin(), clusters_.end(), [](auto* c) { return c->hostWork(); })) {
        logWarning() << "Some clusters compute on the host (see SEISSOL_DEVICE_HOST_SWITCH), or"
                     << "exchange data with clusters that do; no super-timestep gets recorded.";
      }
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

  setUpWatchdog(allClusters);
}

void TimeManager::setUpWatchdog(const std::vector<AbstractTimeCluster*>& clusters) {
  const auto settings = monitoring::Watchdog::settingsFromEnv(seissolInstance_.env());
  // the clusters report and fail after as long without progress as the watchdog
  for (auto* cluster : clusters) {
    cluster->setProgressTimeout(settings.timeout, settings.abort);
  }
  if (!settings.enabled) {
    return;
  }

  watchdog_ = std::make_unique<monitoring::Watchdog>(settings);
  for (auto* cluster : clusters) {
    watchdog_->addProbe("the progress of " + cluster->description(),
                        [cluster] { return cluster->progressValue(); });
    if (const auto* time = cluster->deviceTime(); time != nullptr) {
      // written back by the device
      watchdog_->addProbe("the device time of " + cluster->description(),
                          [time] { return *static_cast<const volatile double*>(time); });
    }
    const bool computes = std::find(clusters_.begin(), clusters_.end(), cluster) != clusters_.end();
    if (computes && parallel::runtime::HostRuntime::scheduling() ==
                        parallel::runtime::HostScheduling::TaskGraph) {
      // in a task graph, the work comes after the actions
      watchdog_->addProbe("the work done of " + cluster->description(),
                          [cluster] { return static_cast<double>(cluster->hostWorkDone()); });
    }
    // each cluster reports on the thread it acts on, if that still runs
    watchdog_->addReport([cluster] { cluster->requestReport(); });
  }
  if (auto* scheduler = haloTransports_->scheduler(); scheduler != nullptr) {
    watchdog_->addCheck([scheduler] { scheduler->checkAsyncErrors(); });
  }
#ifdef SEISSOL_HOST_TASK_GRAPH
  if (auto* engine = haloTransports_->taskEngine(); engine != nullptr) {
    // the main thread may wait for the tasks, and not get to report
    watchdog_->addReport([engine] {
      logWarning(true) << "Watchdog: the halo exchange in the task graph:" << engine->describe();
    });
  }
#endif
  watchdog_->start();
  logInfo() << "A watchdog looks at the progress, and reports"
            << (settings.abort ? "and aborts" : "") << "after" << settings.timeout.count()
            << "s without any (see SEISSOL_WATCHDOG).";
}

void TimeManager::setWatchdogPhase(monitoring::Watchdog::Phase phase, const char* label) {
  if (watchdog_ != nullptr) {
    watchdog_->setPhase(phase, label);
  }
}

void TimeManager::writeClustering() const {
  if (clusteringWriter_.has_value()) {
    clusteringWriter_->write();
  }
}

TimeManager::Updates TimeManager::updates() const {
  Updates updates;
  const auto clusterCount = clusterLayout_.value().globalClusterCount;
  updates.cells.resize(clusterCount);
  updates.duplicateCells.resize(clusterCount);
  updates.dynamicRuptureFaces.resize(clusterCount);
  for (const auto& size : clusterSizes_) {
    const auto steps = static_cast<double>(size.cluster->numberOfTimeSteps());
    updates.cells[size.globalClusterId] += static_cast<double>(size.cells) * steps;
    updates.duplicateCells[size.globalClusterId] +=
        static_cast<double>(size.duplicateCells) * steps;
    updates.dynamicRuptureFaces[size.globalClusterId] +=
        static_cast<double>(size.dynamicRuptureFaces) * steps;
  }
  return updates;
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

  setWatchdogPhase(monitoring::Watchdog::Phase::Stepping, "the time stepping");

  for (auto& cluster : clusters_) {
    cluster->setSyncTime(synchronizationTime);

    // dereference first due to clang-tidy recommendation
    (*cluster).reset();
  }

  communicationManager_->reset(synchronizationTime);

  // clusters that disagree on the number of ticks until the synchronization point would wait for
  // each other forever
  const auto ticks = clusters_.empty() ? 0 : clusters_.front()->getStepsUntilSync();
  const auto agrees = [&](const auto& cluster) { return cluster->getStepsUntilSync() == ticks; };
  const auto& ghostClusters = *communicationManager_->getGhostClusters();
  if (!std::all_of(clusters_.begin(), clusters_.end(), agrees) ||
      !std::all_of(ghostClusters.begin(), ghostClusters.end(), agrees)) {
    logError() << "The clusters disagree on the number of steps until the synchronization point at"
               << synchronizationTime << "s.";
  }

  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());
#ifdef ACL_DEVICE
  const device::DeviceInstance& device = device::DeviceInstance::instance();
  device.api().putProfilingMark("advanceInTime", device::ProfilingColors::Blue);
#endif
  // what the device and the threads do in this epoch
  deviceTimeline_.beginEpoch();
  parallel::runtime::HostRuntime::resetLoopTime();
  const auto epochBegin = std::chrono::steady_clock::now();

  // with a task graph on the host, the actions create the tasks of the epoch; this returns once
  // they are all done
  parallel::runtime::HostRuntime::epoch([this]() {
    actUntilSynchronized();
    // the exchanges of the epoch may still need this thread to get done
    communicationManager_->awaitExchanges();
  });
  if (auto* scheduler = haloTransports_->scheduler(); scheduler != nullptr) {
    // e.g. the exchange engine of actor-mpi stops once it has done the exchanges of the epoch
    scheduler->quiesce();
  }
  // the work on the host that the device may still wait for (e.g. the MPI operations of the
  // actor-mpi transfer mode), and the rest of it
  parallel::runtime::HostActor::instance().drain();
#ifdef ACL_DEVICE
  if (concurrent_) {
    // the clusters have only enqueued their work
    device.api().syncDevice();
  }
  device.api().popLastProfilingMark();
#endif
  // the work on the host for the data of the epoch, before anything reads what it writes (e.g. the
  // samples of the receivers at the synchronization point)
  parallel::runtime::HostActor::instance().drain();
  for (auto& cluster : clusters_) {
    cluster->finishPhase();
  }
  if (recorder_ != nullptr) {
    // the super-timesteps that an epoch has not replayed are not likely to come again (e.g. after
    // the step sizes have changed); their recordings would hold memory of the device to the end
    recorder_->forgetUnused();
    seenSuperStepsBefore_ = std::exchange(seenSuperSteps_, {});
  }
  finishActivity(epochBegin);
  if (monitoring::Trace::recording()) {
    std::ostringstream name;
    name << "epoch to " << synchronizationTime << " s";
    monitoring::Trace::span(monitoring::Trace::name(name.str()),
                            monitoring::Trace::track(this, "time stepping"),
                            monitoring::Trace::time(epochBegin),
                            monitoring::Trace::now());
  }
  setWatchdogPhase(monitoring::Watchdog::Phase::Other, "the end of an epoch");
}

void TimeManager::finishActivity(std::chrono::steady_clock::time_point epochBegin) {
  // the device is done with the epoch; the clusters have added the times of their loops
  for (auto* timer : {superStepTimer_.get(), exchangeTimer_.get()}) {
    if (timer != nullptr) {
      timer->collectAll();
    }
  }
  lastActivity_.device = deviceTimeline_.endEpoch();
  lastActivity_.wall =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - epochBegin).count();
  lastActivity_.threadLoops = parallel::runtime::HostRuntime::loopTime();
  lastActivity_.threads = static_cast<double>(OpenMP::threadCount());

  totalActivity_.device += lastActivity_.device;
  totalActivity_.wall += lastActivity_.wall;
  totalActivity_.threadLoops += lastActivity_.threadLoops;
  totalActivity_.threads = lastActivity_.threads;
}

void TimeManager::actUntilSynchronized() {
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
    throttleHostWork();
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
    watchProgress();
  }
}

void TimeManager::throttleHostWork() {
  // the work on the host whose data the device has written meanwhile
  parallel::runtime::HostActor::instance().progress();
  // with a task graph, the actions do not run too far ahead of their work
  parallel::runtime::HostRuntime::throttle([this]() { communicationManager_->progression(); });
}

void TimeManager::setUpHostScheduling() {
  using parallel::runtime::HostRuntime;
  using parallel::runtime::HostScheduling;
  auto& env = seissolInstance_.env();
  const auto taskSize = hostTaskSize(env);
  if (!useHostTaskGraph(env)) {
    HostRuntime::setScheduling(HostScheduling::Static, taskSize);
    return;
  }
  if (isDeviceOn() && !useHostActor(env)) {
    logWarning() << "A task graph on the host (SEISSOL_HOST_TASK_GRAPH=1) is only available in"
                 << "builds without a device, or with the host actor (SEISSOL_HOST_ACTOR=1);"
                 << "the loops keep a static schedule.";
    HostRuntime::setScheduling(HostScheduling::Static, taskSize);
    return;
  }
  if (!HostRuntime::taskGraphAvailable()) {
    logWarning()
        << "A task graph on the host (SEISSOL_HOST_TASK_GRAPH=1) needs OpenMP 5.1 with "
           "taskloop, masked and detached tasks, which this build does not have (e.g. with "
           "NVHPC, or without OpenMP); the loops keep a static schedule.";
    HostRuntime::setScheduling(HostScheduling::Static, taskSize);
    return;
  }
  if (OpenMP::threadCount() == 1 && Mpi::mpi.size() > 1 &&
      !HostRuntime::detachedTasksInSingleThread()) {
    logWarning() << "With a single thread, the halo exchange in a task graph on the host"
                 << "(SEISSOL_HOST_TASK_GRAPH=1) needs the OpenMP runtime of LLVM 19 or newer;"
                 << "older ones fail an assertion after the first epoch. The loops keep a static"
                 << "schedule.";
    HostRuntime::setScheduling(HostScheduling::Static, taskSize);
    return;
  }
  HostRuntime::setScheduling(HostScheduling::TaskGraph, taskSize);
  HostRuntime::setActionLimit(hostActionLimit(env));
  if (isDeviceOn()) {
    logInfo() << "The work on the host for the device runs in a task graph, while the main thread"
              << "enqueues the work of the device; with tasks of" << taskSize
              << "elements in its loops (SEISSOL_HOST_TASK_SIZE).";
    return;
  }
  logInfo() << "The clusters run their work on the host in a task graph, with tasks of" << taskSize
            << "cells or faces in their loops (SEISSOL_HOST_TASK_SIZE), and at most"
            << HostRuntime::actionLimit() << "actions waiting for their work"
            << "(SEISSOL_HOST_ACTION_LIMIT; 0 for no limit).";
}

void TimeManager::watchProgress() {
  // taking the time costs more than a pass through the loop, and the timeout is in minutes
  constexpr std::size_t PassesPerCheck = 1024;
  if (++waitingPasses_ % PassesPerCheck == 0) {
    // the exchanges the device is past, without waiting
    if (exchangeTimer_ != nullptr) {
      exchangeTimer_->collect();
    }
    for (auto* cluster : clusters_) {
      if (!cluster->synced()) {
        cluster->checkProgress();
      }
    }
  }
}

void TimeManager::completeExchanges(long end) {
  auto* scheduler = haloTransports_->scheduler();
  if (scheduler == nullptr) {
    return;
  }
  // the ghost clusters start the remaining exchanges whose data is complete within the
  // super-timestep
  std::size_t attempts = 0;
  while (!scheduler->launchedBefore(end)) {
    communicationManager_->progression();
    if (++attempts > 100000000) {
      logError() << "The halo exchanges of a super-timestep did not all start.";
    }
  }
}

StepWork TimeManager::takeSteps(const std::vector<PlannedAction>& plan,
                                std::size_t begin,
                                std::size_t end) {
  StepWork work;
  for (auto index = begin; index < end; ++index) {
    const auto& step = plan[index];
    auto* cluster = clusters_[step.cluster];
    throttleHostWork();
    // along the plan, a cluster can only have to wait for the halo exchange
    auto action = cluster->getNextLegalAction();
    while (action == ActorAction::Nothing) {
      communicationManager_->progression();
      parallel::runtime::HostActor::instance().progress();
      watchProgress();
      action = cluster->getNextLegalAction();
    }
    if (action != step.action) {
      logError() << "The cluster" << cluster->identifier() << "is not ready for step" << step.step
                 << "of the time stepping plan.";
    }
    cluster->act();
    work.outputs = work.outputs || cluster->lastStepWork().outputs;
    work.hostWork = work.hostWork || cluster->lastStepWork().hostWork;
  }
  return work;
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
  for (const auto& cluster : planned) {
    largestRate = std::max(largestRate, cluster.timeStepRate);
  }
  const auto plan = planTimeSteps(planned);

  // what a recording of a super-timestep comprises: the clusters, ghost clusters included, and the
  // stream runtimes that their device work goes through (those of the clusters, of their outputs,
  // and of the halo exchange)
  auto* scheduler = haloTransports_->scheduler();
  std::vector<AbstractTimeCluster*> ghosts;
  for (auto& ghost : *communicationManager_->getGhostClusters()) {
    ghosts.push_back(ghost.get());
  }
  std::vector<AbstractTimeCluster*> participants(clusters_.begin(), clusters_.end());
  participants.insert(participants.end(), ghosts.begin(), ghosts.end());
  auto runtimes = scheduler != nullptr ? scheduler->runtimes() : SuperStepRecorder::Runtimes{};
  for (auto* cluster : clusters_) {
    for (auto* runtime : cluster->runtimes()) {
      if (std::find(runtimes.begin(), runtimes.end(), runtime) == runtimes.end()) {
        runtimes.push_back(runtime);
      }
    }
  }

  for (const auto& superTimestep : superTimesteps(planned, plan)) {
    const auto begin = superTimestep.begin;
    const auto end = superTimestep.end;
    const auto superStep = superTimestep.index;
    const auto& steps = superTimestep.steps;

    // decide on the super-timestep before any of its host parts has run
    bool outputsAhead = false;
    bool hostWork = false;
    // a replay repeats the step sizes of its recording, which the key stands for
    bool regularSteps = true;
    SuperStepRecorder::Key key;
    for (std::size_t cluster = 0; cluster < clusters_.size(); ++cluster) {
      outputsAhead = outputsAhead || clusters_[cluster]->outputsAhead(steps[cluster]);
      hostWork = hostWork || clusters_[cluster]->hostWork();
      regularSteps = regularSteps && clusters_[cluster]->regularStepsAhead(steps[cluster]);
      key.push_back(clusters_[cluster]->nextTimeStepSize());
    }
    for (auto* ghost : ghosts) {
      key.push_back(ghost->nextTimeStepSize());
    }
    const bool full = superTimestep.beforeSync;
    const bool replayable = replay_ && full && regularSteps && !outputsAhead && !hostWork;
    const bool reuse = graphMode_ != SuperStepGraphMode::Rebuild;
    const bool replayNow = reuse && replayable && recorder_->has(key);
    // recorded to be replayed later, or else (with the experimental modes) launched once
    const bool seen = seenSuperSteps_.count(key) > 0 || seenSuperStepsBefore_.count(key) > 0;
    const bool keep =
        reuse && replayable && !replayNow && (seen || graphMode_ == SuperStepGraphMode::Always);
    const bool once =
        !replayNow && !keep && replay_ && graphMode_ != SuperStepGraphMode::Replay && !hostWork;

    StepWork work;
    const auto hostStart = std::chrono::steady_clock::now();
    const auto seconds = [](auto from) {
      return std::chrono::duration<double>(std::chrono::steady_clock::now() - from).count();
    };
    const auto superStepEnd = (superStep + 1) * largestRate;
    if (replay_ && scheduler != nullptr) {
      // each super-timestep launches exactly the exchanges whose data it completes, so that
      // recordings and replays comprise the same ones
      scheduler->setHorizon(superStepEnd);
    }
    if (replayNow) {
      recorder_->beginReplay(runtimes);
      for (auto* cluster : clusters_) {
        cluster->setDeviceWork(false);
      }
      if (scheduler != nullptr) {
        scheduler->setLaunching(false);
      }
      work = takeSteps(plan, begin, end);
      completeExchanges(superStepEnd);
      if (scheduler != nullptr) {
        scheduler->setLaunching(true);
      }
      for (auto* cluster : clusters_) {
        cluster->setDeviceWork(true);
      }
      recorder_->replay(key, participants, runtimes, superStepDeviceTimes_.data());
      if (scheduler != nullptr) {
        scheduler->setLatestEvent(recorder_->lastEvent());
      }
      ++replayedSuperSteps_;
      superStepHostTimes_[1] += seconds(hostStart);
    } else if (keep || once) {
      recorder_->beginRecording(participants, runtimes);
      if (scheduler != nullptr) {
        // the recording must not wait for anything from before it: neither for the pending
        // operations, nor for the latest group, which its ghost cluster may acknowledge only inside
        // the recording (e.g. one launched by setHorizon() above); its launch does
        scheduler->restartAfter(ActorEvent());
      }
      work = takeSteps(plan, begin, end);
      completeExchanges(superStepEnd);
      const auto ending = std::chrono::steady_clock::now();
      if (keep) {
        recorder_->endRecording(key);
        recorder_->replay(key, participants, runtimes, &superStepDeviceTimes_[1]);
        ++recordedSuperSteps_;
      } else {
        recorder_->endRecordingAndLaunch(participants, runtimes, &superStepDeviceTimes_[2]);
        ++launchedOnceSuperSteps_;
      }
      superStepHostTimes_[3] += seconds(ending);
      if (scheduler != nullptr) {
        // nothing after the recording may wait for an event recorded inside of it; the streams of
        // the scheduler wait for the replay instead
        scheduler->restartAfter(recorder_->lastEvent());
      }
      superStepHostTimes_[2] += seconds(hostStart);
    } else {
      work = takeSteps(plan, begin, end);
      if (replay_) {
        completeExchanges(superStepEnd);
      }
      if (replayable) {
        seenSuperSteps_.insert(key);
      }
      superStepHostTimes_[0] += seconds(hostStart);
    }

    // the decision in advance has to match what the host parts have decided
    if (work.outputs != outputsAhead) {
      ++mispredictedSuperSteps_;
    }

    ++superSteps_;
    // the replays the device is past, without waiting
    if (superStepTimer_ != nullptr) {
      superStepTimer_->collect();
    }
    if (!full) {
      ++shortenedSuperSteps_;
    } else {
      outputFreeSuperSteps_ += work.outputs ? 0 : 1;
      regularSuperSteps_ += work.irregular() ? 0 : 1;
    }
  }

  // the runtimes may still defer their waits for the latest replay to their next work; they
  // enqueue them now, so that nobody holds an event of the recorder beyond the epoch (the outputs
  // may let go of their runtimes after the last one)
  for (auto* runtime : runtimes) {
    runtime->enqueueDeferredWait();
  }
}

json::Value TimeManager::printComputationTime() {
  actorStateStatisticsManager_.finish();
  auto result = loopStatistics_.printSummary(Mpi::mpi.comm());

  // the loops replayed from recordings of super-timesteps count in the FLOP of the loops, but their
  // time counts as that of the replays
  std::size_t replayed = replayedSuperSteps_ + recordedSuperSteps_ + launchedOnceSuperSteps_;
  MPI_Allreduce(
      MPI_IN_PLACE, &replayed, 1, Mpi::castToMpiType<std::size_t>(), MPI_SUM, Mpi::mpi.comm());
  if (replayed > 0) {
    logInfo() << "The throughput of the loops is not known, as some of their work ran in replays"
              << "of super-timesteps.";
    result["activity"] = printActivity();
    return result;
  }

  // the work of each loop, by the metrics of the flop counter, over the time it took
  const std::vector<std::pair<std::string, std::vector<std::string>>> work{
      {"computeLocalIntegration", {"local"}},
      {"computeNeighboringIntegration",
       {"neighbor", "neighbor-dr", "plasticity-check", "plasticity-yield"}},
      {"computeDynamicRupture", {"dr-frictionlaw-copy", "dr-frictionlaw-interior"}}};
  for (unsigned region = 0; region < loopStatistics_.regionCount(); ++region) {
    const auto& name = loopStatistics_.regionName(region);
    const auto metrics = std::find_if(
        work.begin(), work.end(), [&](const auto& entry) { return entry.first == name; });
    if (metrics == work.end()) {
      continue;
    }
    const auto estimate = seissolInstance_.flopCounter().estimate(metrics->second);
    std::array<double, 3> sums{static_cast<double>(estimate.hardwareFlop),
                               static_cast<double>(estimate.kernelBytes),
                               loopStatistics_.regionTime(region)};
    MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM, Mpi::mpi.comm());
    if (sums[2] <= 0) {
      continue;
    }
    const double flopPerSecond = sums[0] / sums[2];
    const double bytesPerSecond = sums[1] / sums[2];
    logInfo() << name
              << "(throughput while it runs):" << UnitFlopPerS.formatPrefix(flopPerSecond).c_str()
              << "hardware," << UnitBytePerS.formatPrefix(bytesPerSecond).c_str()
              << "kernel bytes (per rank)";
    auto& throughput = result["regions"][name]["throughput"];
    throughput["hardwareFlopPerSecond"] = flopPerSecond;
    throughput["kernelBytesPerSecond"] = bytesPerSecond;
  }

  result["activity"] = printActivity();
  return result;
}

json::Value TimeManager::printActivity() const {
  const auto& total = totalActivity_;
  json::Value activity;

  // the time the threads spent in the loops, of all the time they had
  const double capacity = total.threads * total.wall;
  const auto threadLoops =
      statistics::parallelSummary(capacity > 0 ? total.threadLoops / capacity : 0.0);
  // with a device, the threads may run no loops at all
  if (!DeviceTimeline::active() || threadLoops.max > 0) {
    logInfo() << "Time the threads spent in loops:" << 100 * threadLoops.mean
              << "% (mean over all ranks; minimum" << 100 * threadLoops.min << "%)";
  }
  activity["threadLoops"] = monitoring::toJson(threadLoops);

  if (DeviceTimeline::active()) {
    const auto share = [&](double time) {
      return statistics::parallelSummary(total.device.span > 0 ? time / total.device.span : 0.0);
    };
    const auto computing = share(total.device.computing);
    const auto loops = share(total.device.busyWith(DeviceTimeline::Kind::Loops));
    const auto superSteps = share(total.device.busyWith(DeviceTimeline::Kind::SuperSteps));
    const auto exchanges = share(total.device.busyWith(DeviceTimeline::Kind::Exchanges));
    const auto idle = share(total.device.span - total.device.busy);
    logInfo() << "Time the device computed:" << 100 * computing.mean
              << "% (mean over all ranks; minimum" << 100 * computing.min << "%), in loops"
              << 100 * loops.mean << "% and in replays of super-timesteps" << 100 * superSteps.mean
              << "%; it exchanged the halo" << 100 * exchanges.mean << "% of the time, and was idle"
              << 100 * idle.mean << "% (maximum" << 100 * idle.max << "%)";
    activity["deviceComputing"] = monitoring::toJson(computing);
    activity["deviceLoops"] = monitoring::toJson(loops);
    activity["deviceSuperSteps"] = monitoring::toJson(superSteps);
    activity["deviceExchanges"] = monitoring::toJson(exchanges);
    activity["deviceIdle"] = monitoring::toJson(idle);
  }
  return activity;
}

double TimeManager::getTimeTolerance() const {
  // the same tolerance as the clusters count their steps with
  return ClusterTimes::TickTolerance * clusterLayout_.value().minimumTimestep;
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
  // the watchdog looks at the clusters
  watchdog_.reset();
  // the snapshots copy on a stream of their own
  hostSnapshot_.reset();

  // Nothing may refer to an event of a pool that goes away. The halo exchange holds events of the
  // clusters and of the recorder; the ghost clusters hold events of the halo exchange, of the copy
  // layers and of the recorder; all clusters hold events of the recorder. So they let go of them in
  // this order, before the clusters, the recorder and the halo exchange dispose of their pools.
  auto* scheduler = haloTransports_->scheduler();
  // the timers wait for their events before they give them back
  if (scheduler != nullptr) {
    scheduler->setTimer(nullptr, nullptr);
  }
  if (recorder_ != nullptr) {
    recorder_->setTimer(nullptr);
  }
  exchangeTimer_.reset();
  superStepTimer_.reset();
  deviceTimeline_.release();
  if (scheduler != nullptr) {
    scheduler->forgetEvents();
  }

  // every message sent has to be received somewhere
  std::array<std::size_t, 2> messages{0, 0};
  for (auto& cluster : *communicationManager_->getGhostClusters()) {
    messages[0] += cluster->sentMessages();
    messages[1] += cluster->receivedMessages();
    cluster->finalize();
  }
  // the MPI operations of the actor-mpi transfer mode are done with the transports
  parallel::runtime::HostActor::instance().stopProgressThread();
  for (auto& cluster : clusters_) {
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
  if (scheduler != nullptr && DeviceTimeline::active()) {
    // on the streams of the exchange, from where each group could start, the waits for the
    // other processes included; not within recordings of super-timesteps
    double exchangeTime = exchangeDeviceTime_;
    MPI_Allreduce(MPI_IN_PLACE, &exchangeTime, 1, MPI_DOUBLE, MPI_MAX, Mpi::mpi.comm());
    logInfo() << "Device time of the halo exchange:" << exchangeTime
              << "s (maximum over all ranks)";
  }

  if (scheduler != nullptr) {
    // the events of the halo exchange get reused once nobody waits for them any more
    const auto usage = scheduler->eventUsage();
    std::array<std::size_t, 2> events{usage.pool, usage.peak};
    MPI_Allreduce(MPI_IN_PLACE,
                  events.data(),
                  events.size(),
                  Mpi::castToMpiType<std::size_t>(),
                  MPI_MAX,
                  Mpi::mpi.comm());
    if (events[0] > 0) {
      logInfo() << "Halo exchange events: at most" << events[1] << "of a pool of" << events[0]
                << "referenced at once (maximum over all ranks)";
    }
  }

  if (followPlan_) {
    std::array<std::size_t, 8> superSteps{superSteps_,
                                          shortenedSuperSteps_,
                                          outputFreeSuperSteps_,
                                          regularSuperSteps_,
                                          recordedSuperSteps_,
                                          replayedSuperSteps_,
                                          mispredictedSuperSteps_,
                                          launchedOnceSuperSteps_};
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
              << "; recorded:" << superSteps[4] << "; replayed:" << superSteps[5]
              << "; recorded and launched once:" << superSteps[7] << "(summed over all ranks)";
    if (replay_) {
      auto times = superStepHostTimes_;
      MPI_Allreduce(MPI_IN_PLACE, times.data(), times.size(), MPI_DOUBLE, MPI_MAX, Mpi::mpi.comm());
      logInfo() << "Host time of the super-timesteps: run directly" << times[0] << "s, replayed"
                << times[1] << "s, recorded" << times[2] << "s, of that ending the recordings"
                << times[3] << "s (maximum over all ranks)";
      if (DeviceTimeline::active()) {
        // from where each could start on the device to where it had run; the steps run directly
        // count in the loop statistics
        auto deviceTimes = superStepDeviceTimes_;
        MPI_Allreduce(MPI_IN_PLACE,
                      deviceTimes.data(),
                      deviceTimes.size(),
                      MPI_DOUBLE,
                      MPI_MAX,
                      Mpi::mpi.comm());
        logInfo() << "Device time of the super-timesteps: replayed" << deviceTimes[0]
                  << "s, recorded" << deviceTimes[1] << "s, recorded and launched once"
                  << deviceTimes[2] << "s (maximum over all ranks)";
      }
    }
    if (superSteps[6] > 0) {
      logError() << "For" << superSteps[6]
                 << "super-timesteps, the output samples were not decided correctly in advance.";
    }
  }
  if (messages[0] != messages[1]) {
    logWarning() << "The halo exchange sent and received a different number of messages.";
  }

  // The recordings go before the halo exchange: a graph that has captured operations of its
  // library keeps them in use, e.g. ncclCommDestroy waits until such graphs are destroyed.
  if (recorder_ != nullptr) {
    recorder_->dispose();
  }
  communicationManager_.reset(nullptr);
  haloTransports_.reset(nullptr);
}

bool TimeManager::clustersOnDevice() const {
  return std::all_of(clusters_.begin(), clusters_.end(), [](const auto* cluster) {
    return cluster->getExecutor() == Executor::Device;
  });
}

void TimeManager::synchronizeToHost(const initializer::HostData& data,
                                    [[maybe_unused]] bool fromSnapshot) {
#ifdef ACL_DEVICE
  bool sameExecutor = true;
  for (auto& cluster : clusters_) {
    sameExecutor &= clusters_.front()->getExecutor() == cluster->getExecutor();
  }
  if (!sameExecutor) {
    synchronizeTo(initializer::AllocationPlace::Host);
    return;
  }
  if (fromSnapshot && !data.isEverything() && !data.empty()) {
    if (hostSnapshot_ == nullptr) {
      hostSnapshot_ = std::make_unique<initializer::HostSnapshot>();
    }
    logInfo() << "Synchronizing data... (device->host, through snapshots;" << data.size()
              << "variables)";
    hostSnapshot_->take(data.containers(),
                        device::DeviceInstance::instance().api().getDefaultStream());
    return;
  }
  seissolInstance_.memoryManager().synchronizeToHost(data);
#endif
}

void TimeManager::copySnapshotToHost() {
  if (hostSnapshot_ != nullptr) {
    hostSnapshot_->copyToHost();
  }
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
    auto* stream = device::DeviceInstance::instance().api().getDefaultStream();
    for (auto& cluster : clusters_) {
      cluster->synchronizeTo(place, stream);
    }
    device::DeviceInstance::instance().api().syncDefaultStreamWithHost();
  }
#endif
}

} // namespace seissol::solver

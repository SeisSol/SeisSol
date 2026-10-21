// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)
// SPDX-FileContributor: Sebastian Rettenberger

#include "CellCluster.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "Common/Marker.h"
#include "Config.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Receiver.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/DeviceTimeline.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "SeisSol.h"
#include "Solver/Settings.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Actor/StepParams.h"

#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <string>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/DeviceGraph.h"

#include <Device/AbstractAPI.h>
#endif

namespace seissol::solver {

template <typename Cfg>
CellCluster<Cfg>::CellCluster(unsigned int clusterId,
                              unsigned int globalClusterId,
                              unsigned int profilingId,
                              const SimulationSettings& settings,
                              HaloType layerType,
                              double maxTimeStepSize,
                              long timeStepRate,
                              bool printProgress,
                              CompoundGlobalData<Cfg> globalData,
                              LTS::Layer* clusterData,
                              seissol::SeisSol& seissolInstance,
                              LoopStatistics* loopStatistics,
                              DeviceTimeline* deviceTimeline,
                              ActorStateStatistics* actorStateStatistics)
    : CellClusterInterface(
          maxTimeStepSize, timeStepRate, seissolInstance.executionPlace(clusterData->size())),
      // cluster ids
      settings_(settings), seissolInstance_(seissolInstance),
      configBoundary_(seissolInstance.parameters().model.configs()), streamRuntime_(4),
      globalData_(globalData), clusterData_(clusterData),
      sourceCluster_(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      loopStatistics_(loopStatistics),
      deviceLoops_(loopStatistics, deviceTimeline, DeviceTimeline::Kind::Loops),
      actorStateStatistics_(actorStateStatistics),
      conditionalCounterHost_(1, seissol::memory::Memkind::Standard),
      conditionalCounterDevice_(1,
                                isDeviceOn() ? seissol::memory::Memkind::DeviceGlobalMemory
                                             : seissol::memory::Memkind::Standard),
      layerType_(layerType), printProgress_(printProgress), clusterId_(clusterId),
      globalClusterId_(globalClusterId), profilingId_(profilingId) {
  // assert all pointers are valid
  assert(clusterData_ != nullptr);
  assert(globalData_.onHost != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalData_.onDevice != nullptr);
  }

  // set timings to zero
  receiverTime_ = 0;

  spacetimeKernel_.setGlobalData(globalData);
  timeKernel_.setGlobalData(globalData);
  localKernel_.setGlobalData(globalData);
  localKernel_.setInitConds(&seissolInstance_.memoryManager().initialConditions(configIdOf<Cfg>()));
  localKernel_.setGravitationalAcceleration(seissolInstance_.gravitationSetup().acceleration);
  neighborKernel_.setGlobalData(globalData);
  if constexpr (seissol::isDeviceOn()) {
    // the device kernels of the faces between configurations read the constants of both
    for (const auto config : seissolInstance_.parameters().model.configs()) {
      dispatchConfig(config, [&](auto configCfg) {
        using ConfigCfg = decltype(configCfg);
        configBoundary_.template setDeviceGlobalData<ConfigCfg>(
            seissolInstance_.memoryManager().template globalData<ConfigCfg>().onDevice);
      });
    }
  }

  computeFlops();

  regionComputeLocalIntegration_ = loopStatistics_->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration_ =
      loopStatistics_->getRegion("computeNeighboringIntegration");
  regionComputePointSources_ = loopStatistics_->getRegion("computePointSources");

  conditionalCounterHost_[0] = 0;
  conditionalCounterDevice_.copyFrom(conditionalCounterHost_);

  perfHandle_[static_cast<std::size_t>(ComputePart::Local)] =
      seissolInstance.flopCounter().addMetric("local", "WP");
  perfHandle_[static_cast<std::size_t>(ComputePart::Neighbor)] =
      seissolInstance.flopCounter().addMetric("neighbor", "WP");
  perfHandle_[static_cast<std::size_t>(ComputePart::DRNeighbor)] =
      seissolInstance.flopCounter().addMetric("neighbor-dr", "DR");
  perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityCheck)] =
      seissolInstance.flopCounter().addMetric("plasticity-check", "PL");
  perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityYield)] =
      seissolInstance.flopCounter().addMetric("plasticity-yield", "PL");

  const auto* cellInfo = clusterData_->var<LTS::CellInformation>();
  for (std::size_t i = 0; i < clusterData_->size(); ++i) {
    if (cellInfo[i].plasticityEnabled) {
      ++numPlasticCells_;
    }
  }
}

template <typename Cfg>
void CellCluster<Cfg>::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster_ = std::move(sourceCluster);
  if (sourceCluster_.device != nullptr) {
    sourceCluster_.device->setClock(clock_.device());
  }
}

template <typename Cfg>
void CellCluster<Cfg>::writeReceivers(const StepParams& params,
                                      const kernels::ReceiverCluster::Sampling& sampling) {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster_ != nullptr) {
    if (runTimeOutputs_) {
      receiverCluster_->sampleAtRunTime(
          clock_.device(), params.time, params.timeStepSize, executor_, streamRuntime_);
    } else {
      receiverCluster_->sample(
          sampling, params.time, params.timeStepSize, executor_, streamRuntime_);
    }
  }
}

template <typename Cfg>
void CellCluster<Cfg>::computeSources(const StepParams& params) {
#ifdef ACL_DEVICE
  device_.api().putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
  SCOREP_USER_REGION("computeSources", SCOREP_USER_REGION_TYPE_FUNCTION)

  // Return when point sources not initialized. This might happen if there
  // are no point sources on this rank.
  auto* pointSourceCluster = [&]() -> kernels::PointSourceCluster* {
    if (executor_ == Executor::Device) {
      return sourceCluster_.device.get();
    } else {
      return sourceCluster_.host.get();
    }
  }();

  if (pointSourceCluster != nullptr) {
    if (executor_ == Executor::Device) {
      const auto loop = deviceLoops_.begin(streamRuntime_.stream());
      pointSourceCluster->addTimeIntegratedPointSources(
          params.time, params.timeStepSize, streamRuntime_);
      deviceLoops_.end(loop,
                       streamRuntime_.stream(),
                       regionComputePointSources_,
                       pointSourceCluster->size(),
                       profilingId_);
    } else {
      const auto loopBegin = LoopStatistics::begin();
      pointSourceCluster->addTimeIntegratedPointSources(
          params.time, params.timeStepSize, streamRuntime_);
      loopStatistics_->end(
          regionComputePointSources_, pointSourceCluster->size(), profilingId_, loopBegin);
    }
  }
#ifdef ACL_DEVICE
  device_.api().popLastProfilingMark();
#endif
}

template <typename Cfg>
void CellCluster<Cfg>::computeLocalIntegration(const StepParams& params) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  const auto loopBegin = LoopStatistics::begin();

  real* const* stepIntegrals = clusterData_->var<LTS::StepIntegrals>(Cfg());
  real* const* accumulatedIntegrals = clusterData_->var<LTS::AccumulatedIntegrals>(Cfg());
  real* const* derivatives = clusterData_->var<LTS::Derivatives>(Cfg());

  const kernels::LocalTmp<Cfg> tmp(seissolInstance_.gravitationSetup().acceleration);

  const auto timeStepWidth = params.timeStepSize;
  const auto startTime = params.time;
  const auto resetBuffers = params.resetBuffers;
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

  parallel::runtime::HostRuntime::forRange(
      clusterData_->size(), [&](std::size_t begin, std::size_t end) {
        // local integration buffer
        alignas(Alignment) real integrationBuffer[kernels::SolverOf<Cfg>::IntegralsSize];

        auto rangeTmp = tmp;

        for (std::size_t cell = begin; cell < end; ++cell) {
          auto data = clusterData_->cellRef<Cfg>(cell);

          // pointer for the call of the ADER-function
          real* bufferPointer = nullptr;
          if (data.template get<LTS::CellInformation>().ltsSetup.hasBuffer(
                  BufferType::StepIntegrals)) {
            // assert presence of the buffer
            assert(stepIntegrals[cell] != nullptr);

            bufferPointer = stepIntegrals[cell];
          } else {
            // work on local buffer
            bufferPointer = integrationBuffer;
          }

          spacetimeKernel_.computeAder(integrationCoeffs.data(),
                                       timeStepWidth,
                                       data,
                                       rangeTmp,
                                       bufferPointer,
                                       derivatives[cell],
                                       true);

          // Compute local integrals (including local boundary conditions)
          localKernel_.computeIntegral(bufferPointer, data, rangeTmp, startTime, timeStepWidth);

          for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
            auto& curFaceDisplacements = data.template get<LTS::FaceDisplacements>()[face];
            // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
            if (curFaceDisplacements != nullptr &&
                data.template get<LTS::CellInformation>().faceTypes[face] !=
                    FaceType::FreeSurfaceGravity) {
              kernel::addVelocity<Cfg> addVelocityKrnl;

              addVelocityKrnl.bindGlobals(*globalData_.onHost);
              addVelocityKrnl.faceDisplacement = data.template get<LTS::FaceDisplacements>()[face];
              addVelocityKrnl.I = bufferPointer;
              addVelocityKrnl.execute(face);
            }
          }

          // We've used a step integral so far -> accumulate update if needed.
          if (data.template get<LTS::CellInformation>().ltsSetup.hasBuffer(
                  BufferType::AccumulatedIntegrals)) {
            assert(accumulatedIntegrals[cell] != nullptr);

            if (resetBuffers) {
              std::memcpy(accumulatedIntegrals[cell],
                          bufferPointer,
                          kernels::SolverOf<Cfg>::IntegralsSize * sizeof(real));
            } else {
#pragma omp simd
              for (std::size_t dof = 0; dof < kernels::SolverOf<Cfg>::IntegralsSize; ++dof) {
                accumulatedIntegrals[cell][dof] += bufferPointer[dof];
              }
            }
          }
        }
      });

  loopStatistics_->end(
      regionComputeLocalIntegration_, clusterData_->size(), profilingId_, loopBegin);
}

template <typename Cfg>
void CellCluster<Cfg>::computeLocalIntegrationDevice(SEISSOL_GPU_PARAM const StepParams& params) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device_.api().putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  const auto loop = deviceLoops_.begin(streamRuntime_.stream());

  auto& dataTable = clusterData_->getConditionalTable<inner_keys::Wp>();
  auto& indicesTable = clusterData_->getConditionalTable<inner_keys::Indices>();

  kernels::LocalTmp<Cfg> tmp(seissolInstance_.gravitationSetup().acceleration);

  const double timeStepWidth = params.timeStepSize;
  const bool resetBuffers = params.resetBuffers;
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

  const ComputeGraphType graphType =
      resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, true);
  streamRuntime_.runGraph(
      computeGraphKey,
      *clusterData_,
      [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        spacetimeKernel_.computeBatchedAder(integrationCoeffs.data(),
                                            timeStepWidth,
                                            *clusterData_,
                                            tmp,
                                            dataTable,
                                            true,
                                            streamRuntime);

        localKernel_.computeBatchedIntegral(dataTable, indicesTable, timeStepWidth, streamRuntime);

        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          const ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
          if (dataTable.find(key) != dataTable.end()) {
            auto& entry = dataTable[key];
            // NOTE: the integrated velocities are not stored separately; the recorded pointers
            // point into the integrated dofs, at the first velocity column
            // (model::MaterialT::VelocityOffset).

            kernel::gpu_addVelocity<Cfg> displacementKrnl;
            displacementKrnl.faceDisplacement =
                entry.get<real*>(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
            displacementKrnl.integratedVelocities = const_cast<const real**>(
                entry.get<real*>(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
            displacementKrnl.bindGlobals(*globalData_.onDevice);

            // Note: this kernel doesn't require tmp. memory
            displacementKrnl.numElements =
                entry.get<real*>(inner_keys::Wp::Id::FaceDisplacement)->getSize();
            displacementKrnl.streamPtr = streamRuntime_.stream();
            displacementKrnl.execute(face);
          }
        }

        const ConditionalKey key =
            ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
        if (dataTable.find(key) != dataTable.end()) {
          auto& entry = dataTable[key];

          if (resetBuffers) {
            device_.algorithms().streamBatchedData(
                const_cast<const real**>(
                    (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get<real*>(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I<Cfg>::Size,
                (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime_.stream());
          } else {
            device_.algorithms().accumulateBatchedData(
                const_cast<const real**>(
                    (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get<real*>(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I<Cfg>::Size,
                (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime_.stream());
          }
        }
      },
      isRecurringTimestep(timeStepWidth));

  // depends on the current time, and therefore cannot be replayed from the graph above; it neither
  // reads nor writes what the graph computes after the local integral
  localKernel_.evaluateBatchedTimeDependentBc(
      dataTable, indicesTable, *clusterData_, clock_.host(), timeStepWidth, streamRuntime_);

  deviceLoops_.end(loop,
                   streamRuntime_.stream(),
                   regionComputeLocalIntegration_,
                   clusterData_->size(),
                   profilingId_);
  device_.api().popLastProfilingMark();
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

template <typename Cfg>
void CellCluster<Cfg>::computeNeighboringIntegration(const StepParams& params) {
  if (settings_.integrate) {
    if (settings_.plasticity) {
      computeNeighboringIntegrationImplementation<true, true>(params);
    } else {
      computeNeighboringIntegrationImplementation<false, true>(params);
    }
  } else {
    if (settings_.plasticity) {
      computeNeighboringIntegrationImplementation<true, false>(params);
    } else {
      computeNeighboringIntegrationImplementation<false, false>(params);
    }
  }
}

template <typename Cfg>
void CellCluster<Cfg>::computeNeighboringIntegrationDevice(
    SEISSOL_GPU_PARAM const StepParams& params) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  device_.api().putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  const auto loop = deviceLoops_.begin(streamRuntime_.stream());

  const double timeStepWidth = params.timeStepSize;
  auto& table = clusterData_->getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
  const auto subtimeCoeffs = timeBasis.integrate(
      params.subTimeStart, timeStepWidth + params.subTimeStart, params.neighborTimeStepSize);

  seissol::kernels::TimeCommon<Cfg>::computeBatchedIntegrals(
      timeKernel_, timeCoeffs.data(), subtimeCoeffs.data(), table, streamRuntime_);
  if (!configBoundary_.empty()) {
    configBoundary_.setIntervals(timeStepWidth, params.subTimeStart, params.neighborTimeStepSize);
    configBoundary_.computeBatchedIntegrals(table, streamRuntime_);
  }

  const ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime_.runGraph(
      computeGraphKey,
      *clusterData_,
      [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        neighborKernel_.computeBatchedNeighborsIntegral(table, streamRuntime);
      },
      true);

  if (settings_.plasticity) {
    auto plasticityGraphKey = initializer::GraphKey(ComputeGraphType::Plasticity, timeStepWidth);
    auto* plasticity =
        clusterData_->var<LTS::Plasticity>(Cfg(), seissol::initializer::AllocationPlace::Device);
    auto* isAdjustableVector =
        clusterData_->var<LTS::FlagScratch>(seissol::initializer::AllocationPlace::Device);
    auto& indicesTable = clusterData_->getConditionalTable<inner_keys::Indices>();
    streamRuntime_.runGraph(
        plasticityGraphKey,
        *clusterData_,
        [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
          seissol::kernels::Plasticity<Cfg>::computePlasticityBatched(
              timeStepWidth,
              seissolInstance_.parameters().model.tv,
              globalData_.onDevice,
              table,
              indicesTable,
              plasticity,
              conditionalCounterDevice_.data(),
              isAdjustableVector,
              streamRuntime);
        },
        isRecurringTimestep(timeStepWidth));
  }

  if (settings_.integrate) {
    const ConditionalKey key = ConditionalKey(*KernelNames::Time);
    if (table.find(key) != table.end()) {
      auto entry = table.at(key);
      device_.algorithms().accumulateBatchedData(
          const_cast<const real**>(
              (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr()),
          (entry.get<real*>(inner_keys::Wp::Id::Integrals))->getDeviceDataPtr(),
          tensor::Q<Cfg>::Size,
          (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize(),
          streamRuntime_.stream());
    }
  }

  device_.api().popLastProfilingMark();
  deviceLoops_.end(loop,
                   streamRuntime_.stream(),
                   regionComputeNeighboringIntegration_,
                   clusterData_->size(),
                   profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

template <typename Cfg>
void CellCluster<Cfg>::computeLocalIntegrationFlops() {
  auto& estimate = estimate_[static_cast<int>(ComputePart::Local)];
  estimate = PerformanceEstimate{};

  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    estimate += spacetimeKernel_.metrics();
    estimate += localKernel_.metrics(cellInformation[cell].faceTypes);

    // Contribution from displacement/integrated displacement
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::FreeSurfaceGravity) {
        estimate += GravitationalFreeSurfaceBc<Cfg>::metrics(face);
      }
    }
  }
}

template <typename Cfg>
void CellCluster<Cfg>::computeNeighborIntegrationFlops() {
  auto& estimateRegular = estimate_[static_cast<int>(ComputePart::Neighbor)];
  auto& estimateDR = estimate_[static_cast<int>(ComputePart::DRNeighbor)];

  estimateRegular = PerformanceEstimate{};
  estimateDR = PerformanceEstimate{};

  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* drMapping = clusterData_->var<LTS::DRMapping>(Cfg());
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    const auto [cellRegular, cellDR] = neighborKernel_.metrics(
        cellInformation[cell].faceTypes, cellInformation[cell].faceRelations, drMapping[cell]);

    estimateRegular += cellRegular;
    estimateDR += cellDR;
  }
}

template <typename Cfg>
void CellCluster<Cfg>::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();

  const auto [check, yield] = seissol::kernels::Plasticity<Cfg>::metrics();
  estimate_[static_cast<int>(ComputePart::PlasticityCheck)] = check;
  estimate_[static_cast<int>(ComputePart::PlasticityYield)] = yield;
}

template <typename Cfg>
ActResult CellCluster<Cfg>::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

template <typename Cfg>
void CellCluster<Cfg>::handleNeighborPrediction(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
template <typename Cfg>
void CellCluster<Cfg>::handleNeighborCorrection(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
template <typename Cfg>
bool CellCluster<Cfg>::hostWork() const {
  return executor_ == Executor::Host || hasDifferentExecutorNeighbor();
}

template <typename Cfg>
bool CellCluster<Cfg>::outputsAhead(long steps) const {
  if (clusterData_->size() == 0 || receiverCluster_ == nullptr || runTimeOutputs_) {
    return false;
  }
  // the times and receiver times the host parts of the next steps start with
  double time = ct_.correctionTime;
  double receiverTime = receiverTime_;
  for (long step = 0; step < steps; ++step) {
    const double timeStepSize =
        ct_.timeStepSize(time, ct_.stepsSinceLastSync + step * ct_.timeStepRate, syncTime_);
    const auto sampling = receiverCluster_->planSampling(receiverTime, time, timeStepSize);
    if (sampling.due) {
      return true;
    }
    receiverTime = sampling.nextTime;
    time += timeStepSize;
  }
  return false;
}

template <typename Cfg>
StepWork CellCluster<Cfg>::prepare(ActorAction action) {
  const auto params = stepParams();
  StepWork work;
  work.hostWork = hostWork();

  if (action == ActorAction::Predict) {
    // without a device, the clock is up to date on the host; it has to agree with the time
    assert(isDeviceOn() || *clock_.host() == params.time);

    if (clusterData_->size() > 0) {
      if (receiverCluster_ != nullptr && !runTimeOutputs_) {
        receiverSampling_ =
            receiverCluster_->planSampling(receiverTime_, params.time, params.timeStepSize);
        receiverTime_ = receiverSampling_.nextTime;
        work.outputs = receiverSampling_.due;
      }
      incrementPerformanceMetrics(ComputePart::Local);
    }
  }

  if (action == ActorAction::Correct) {
    incrementPerformanceMetrics(ComputePart::Neighbor);
    incrementPerformanceMetrics(ComputePart::DRNeighbor);
    if (settings_.plasticity) {
      // here instead of with the device work, which a replay of a recording skips; the yielding
      // cells get counted on the device, and added at the synchronization points
      seissolInstance_.flopCounter().incrementMetric(
          perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityCheck)],
          estimate_[static_cast<std::size_t>(ComputePart::PlasticityCheck)] * numPlasticCells_);
    }

    if (printProgress_) {
      const auto nextCorrectionSteps = ct_.nextCorrectionSteps();
      if (((nextCorrectionSteps / timeStepRate_) % 100) == 0) {
        logInfo() << "Max cluster / LTS cycle updates since sync: " << nextCorrectionSteps
                  << " at time " << ct_.nextCorrectionTime(syncTime_);
      }
    }
  }
  return work;
}

template <typename Cfg>
void CellCluster<Cfg>::predict() {
  assert(state_ == ActorState::Corrected);
  const auto params = stepParams();
  deviceLoops_.collect();

  if (clusterData_->size() == 0) {
    return;
  }

  if (executor_ == Executor::Device) {
    writeReceivers(params, receiverSampling_);
    computeLocalIntegrationDevice(params);
    computeSources(params);
  } else {
    // the prediction changes what the neighbors read; the work takes the sampling of its step
    // along, as the decisions about the next steps may come before it runs
    hostRuntime_.action({/*readsNeighbors=*/false, /*writesShared=*/true},
                        [this, params, sampling = receiverSampling_]() {
                          writeReceivers(params, sampling);
                          computeLocalIntegration(params);
                          computeSources(params);
                        });
  }

  if (hasDifferentExecutorNeighbor()) {
    auto other = executor_ == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                               : seissol::initializer::AllocationPlace::Device;
    clusterData_->varSynchronizeTo<LTS::Buffers>(other, streamRuntime_.stream());
  }

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

template <typename Cfg>
void CellCluster<Cfg>::correct() {
  assert(state_ == ActorState::Predicted);
  const auto params = stepParams();
  deviceLoops_.collect();

  if (executor_ == Executor::Device) {
    computeNeighboringIntegrationDevice(params);
  } else {
    // the correction reads what the neighbors share, and changes the DOFs only
    hostRuntime_.action({/*readsNeighbors=*/true, /*writesShared=*/false},
                        [this, params]() { computeNeighboringIntegration(params); });
  }

  // the time of the cluster advances by the same step after the correction
  clock_.advance(params.timeStepSize, streamRuntime_);

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

template <typename Cfg>
void CellCluster<Cfg>::timeSet(double time) {
  clock_.set(time, streamRuntime_);
}

template <typename Cfg>
ActorEvent CellCluster<Cfg>::recordActionEvent() {
  return ActorEvent(streamRuntime_.eventRecord());
}

template <typename Cfg>
void CellCluster<Cfg>::waitForEvent(SEISSOL_GPU_PARAM void* event) {
#ifdef ACL_DEVICE
  if (executor_ == Executor::Host) {
    // the host kernels run right away, so the host has to wait
    device_.api().syncEventWithHost(event);
  } else {
    device_.api().syncStreamWithEvent(streamRuntime_.stream(), event);
  }
#endif
}

template <typename Cfg>
void CellCluster<Cfg>::incrementPerformanceMetrics(ComputePart part) {
  seissolInstance_.flopCounter().incrementMetric(perfHandle_[static_cast<std::size_t>(part)],
                                                 estimate_[static_cast<std::size_t>(part)]);
}

template <typename Cfg>
unsigned int CellCluster<Cfg>::getClusterId() const {
  return clusterId_;
}

template <typename Cfg>
std::size_t CellCluster<Cfg>::layerId() const {
  return clusterData_->id();
}

template <typename Cfg>
unsigned int CellCluster<Cfg>::getGlobalClusterId() const {
  return globalClusterId_;
}

template <typename Cfg>
HaloType CellCluster<Cfg>::getLayerType() const {
  return layerType_;
}
template <typename Cfg>
void CellCluster<Cfg>::setTime(double time) {
  AbstractTimeCluster::setTime(time);
  this->receiverTime_ = time;
  if (receiverCluster_ != nullptr) {
    receiverCluster_->setNextSampleTime(time);
  }
}

template <typename Cfg>
void CellCluster<Cfg>::setRunTimeOutputs(bool runTimeOutputs) {
  runTimeOutputs_ = runTimeOutputs;
}

template <typename Cfg>
void* CellCluster<Cfg>::outputStream() {
  return receiverCluster_ != nullptr ? receiverCluster_->sampleStream() : nullptr;
}

template <typename Cfg>
void CellCluster<Cfg>::finalize() {
  sourceCluster_.host.reset(nullptr);
  sourceCluster_.device.reset(nullptr);
  AbstractTimeCluster::finalize();
  deviceLoops_.release();
  clock_.dispose();
  streamRuntime_.dispose();

  logDebug() << "#(time steps):" << numberOfTimeSteps_;
}

template <typename Cfg>
template <bool UsePlasticity, bool IntegrateOutput>
void CellCluster<Cfg>::computeNeighboringIntegrationImplementation(const StepParams& params) {
  const auto clusterSize = clusterData_->size();
  if (clusterSize == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  const auto loopBegin = LoopStatistics::begin();

  const auto* faceNeighbors = clusterData_->var<LTS::FaceNeighbors>();
  const auto* drMapping = clusterData_->var<LTS::DRMapping>(Cfg());
  const auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* plasticity = clusterData_->var<LTS::Plasticity>(Cfg());
  auto* pstrain = clusterData_->var<LTS::PStrain>(Cfg());

  // the cells whose plasticity yields, counted per range of cells and summed up (without
  // plasticity, neither counter changes)
  // NOLINTNEXTLINE(misc-const-correctness)
  std::atomic<std::size_t> numberOfTetsWithPlasticYielding{0};

  const auto tV = seissolInstance_.parameters().model.tv;

  const auto timestep = params.timeStepSize;
  const auto oneMinusIntegratingFactor =
      seissol::kernels::Plasticity<Cfg>::computeRelaxTime(tV, timestep);

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs = timeBasis.integrate(
      params.subTimeStart, timestep + params.subTimeStart, params.neighborTimeStepSize);

  const bool configBoundary = !configBoundary_.empty();
  if (configBoundary) {
    configBoundary_.setIntervals(timestep, params.subTimeStart, params.neighborTimeStepSize);
  }

  parallel::runtime::HostRuntime::forRange(clusterSize, [&](std::size_t begin, std::size_t end) {
    std::array<real*, Cell::NumFaces> timeIntegrated{};
    std::array<real*, Cell::NumFaces> faceNeighborsPrefetch{};
    // NOLINTNEXTLINE(misc-const-correctness)
    std::size_t yielding = 0;

    for (std::size_t cell = begin; cell < end; ++cell) {
      auto data = clusterData_->cellRef<Cfg>(cell);

      // Scratch for the neighbours whose time integral has to be computed here.
      // Written before it is read, so it needs no initialisation; the frame
      // holds it for the whole loop, one copy per range.
      alignas(Alignment)
          real integrationBuffer[Cell::NumFaces][kernels::SolverOf<Cfg>::IntegralsSize];
      std::array<real*, Cell::NumFaces> integrationBuffers{};
      for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
        integrationBuffers[i] = integrationBuffer[i];
      }

      seissol::kernels::TimeCommon<Cfg>::computeIntegrals(timeKernel_,
                                                          data.template get<LTS::CellInformation>(),
                                                          timeCoeffs.data(),
                                                          subtimeCoeffs.data(),
                                                          faceNeighbors[cell],
                                                          integrationBuffers,
                                                          timeIntegrated);
      if (configBoundary) {
        configBoundary_.computeIntegrals(data.template get<LTS::CellInformation>(),
                                         data.template get<LTS::NormalStress>(),
                                         faceNeighbors[cell],
                                         integrationBuffers,
                                         timeIntegrated);
      }

      faceNeighborsPrefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                     ? static_cast<real*>(faceNeighbors[cell][1])
                                     : drMapping[cell][1].godunov;
      faceNeighborsPrefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                     ? static_cast<real*>(faceNeighbors[cell][2])
                                     : drMapping[cell][2].godunov;
      faceNeighborsPrefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                     ? static_cast<real*>(faceNeighbors[cell][3])
                                     : drMapping[cell][3].godunov;

      // fourth face's prefetches: the next cell of the range, which comes next on this thread
      if (cell + 1 < end) {
        faceNeighborsPrefetch[3] =
            (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
                ? static_cast<real*>(faceNeighbors[cell + 1][0])
                : drMapping[cell + 1][0].godunov;
      } else {
        faceNeighborsPrefetch[3] = static_cast<real*>(faceNeighbors[cell][3]);
      }

      neighborKernel_.computeNeighborsIntegral(data, timeIntegrated, faceNeighborsPrefetch);

      if constexpr (UsePlasticity) {
        if (data.template get<LTS::CellInformation>().plasticityEnabled) {
          yielding +=
              seissol::kernels::Plasticity<Cfg>::computePlasticity(oneMinusIntegratingFactor,
                                                                   timestep,
                                                                   tV,
                                                                   globalData_.onHost,
                                                                   &plasticity[cell],
                                                                   data.template get<LTS::Dofs>(),
                                                                   pstrain[cell]);
        }
      }
      if constexpr (IntegrateOutput) {
        auto* __restrict integral = data.template get<LTS::Integrals>();
        const auto* __restrict dofs = data.template get<LTS::Dofs>();

// only first-order time integration for the output here
#pragma omp simd
        for (std::size_t dof = 0; dof < tensor::Q<Cfg>::size(); ++dof) {
          integral[dof] += timestep * dofs[dof];
        }
      }
    }

    if constexpr (UsePlasticity) {
      numberOfTetsWithPlasticYielding.fetch_add(yielding, std::memory_order_relaxed);
    }
  });

  if constexpr (UsePlasticity) {
    conditionalCounterHost_[0] += numberOfTetsWithPlasticYielding.load(std::memory_order_relaxed);
  }

  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterSize, profilingId_, loopBegin);
}

template <typename Cfg>
void CellCluster<Cfg>::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor_ == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor_ == Executor::Host)) {
      clusterData_->synchronizeTo(place, stream);
    }
  }
}

template <typename Cfg>
void CellCluster<Cfg>::finishPhase() {
  // the device is past all loops of the epoch
  deviceLoops_.collectAll();

  const auto cells = conditionalCounterHost_[0];
  seissolInstance_.flopCounter().incrementMetric(
      perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityYield)],
      estimate_[static_cast<std::size_t>(ComputePart::PlasticityYield)] * cells);

  conditionalCounterHost_.copyFrom(conditionalCounterDevice_);
  const auto cells2 = conditionalCounterHost_[0];
  seissolInstance_.flopCounter().incrementMetric(
      perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityYield)],
      estimate_[static_cast<std::size_t>(ComputePart::PlasticityYield)] * cells2);

  conditionalCounterHost_[0] = 0;
  conditionalCounterDevice_.copyFrom(conditionalCounterHost_);
}

template <typename Cfg>
std::string CellCluster<Cfg>::description() const {
  const auto identifier = clusterData_->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  // the time cluster, and the configuration where there can be several
  return "compute-" + haloStr + "-" + std::to_string(identifier.lts) +
         (builtConfigCount() > 1 ? "-config" + std::to_string(identifier.config) : std::string());
}

#define SEISSOL_INSTANTIATE(Cfg) template class CellCluster<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::solver

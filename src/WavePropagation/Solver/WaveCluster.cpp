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

#include "WaveCluster.h"

#include "Alignment.h"
#include "Config.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Constants.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Instrumentation.h"
#include "Core/Common/Marker.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/LtsSetup.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/HostRuntime.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Parallel/Trace.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Compute/CellStage.h"
#include "GeneratedCode/kernel.h"
#include "WavePropagation/Output/Receiver.h"
#include "WavePropagation/Solver/Interface.h"
#include "WavePropagation/Solver/LinearCK/GravitationalFreeSurfaceBC.h"
#include "WavePropagation/Solver/Solver.h"
#include "WavePropagation/Solver/SolverSelector.h"
#include "WavePropagation/Solver/TimeCommon.h"

#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <functional>
#include <memory>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Core/Memory/Batch/ConditionalKey.h"
#include "Core/Memory/Batch/DeviceGraph.h"
#include "Core/Memory/Batch/EncodedConstants.h"
#include "GeneratedCode/tensor.h"

#include <Device/AbstractAPI.h>
#endif

namespace seissol::solver {

template <typename Cfg>
WaveCluster<Cfg>::WaveCluster(unsigned int clusterId,
                              unsigned int globalClusterId,
                              unsigned int profilingId,
                              double maxTimeStepSize,
                              long timeStepRate,
                              Executor executor,
                              bool printProgress,
                              LTS::Layer* clusterData,
                              const Setup& setup,
                              LoopStatistics* loopStatistics,
                              DeviceTimeline* deviceTimeline,
                              ActorStateStatistics* actorStateStatistics)
    : CellCluster(clusterId,
                  globalClusterId,
                  profilingId,
                  maxTimeStepSize,
                  timeStepRate,
                  executor,
                  printProgress,
                  clusterData,
                  loopStatistics,
                  deviceTimeline,
                  actorStateStatistics),
      run_(setup.run), configBoundary_(setup.configs),
      globalData_(setup.globalData.template get<Cfg>()),
      gravitationalAcceleration_(setup.gravitationalAcceleration) {
  assert(globalData_.onHost != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalData_.onDevice != nullptr);
  }

  spacetimeKernel_.setGlobalData(globalData_);
  timeKernel_.setGlobalData(globalData_);
  localKernel_.setGlobalData(globalData_);
  localKernel_.setInitConds(setup.initialConditions);
  localKernel_.setGravitationalAcceleration(gravitationalAcceleration_);
  neighborKernel_.setGlobalData(globalData_);
  if constexpr (seissol::isDeviceOn()) {
    // the device kernels of the faces between configurations read the constants of both
    for (const auto config : setup.configs) {
      dispatchConfig(config, [&](auto configCfg) {
        using ConfigCfg = decltype(configCfg);
        configBoundary_.template setDeviceGlobalData<ConfigCfg>(
            setup.globalData.template get<ConfigCfg>().onDevice);
      });
    }
  }

  computeFlops();

  regionComputeLocalIntegration_ = loopStatistics_->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration_ =
      loopStatistics_->getRegion("computeNeighboringIntegration");
  traceLocalIntegration_ = monitoring::Trace::name("computeLocalIntegration");
  traceNeighboringIntegration_ = monitoring::Trace::name("computeNeighboringIntegration");

  perfHandle_[static_cast<std::size_t>(ComputePart::Local)] =
      run_.flopCounter().addMetric("local", "WP");
  perfHandle_[static_cast<std::size_t>(ComputePart::Neighbor)] =
      run_.flopCounter().addMetric("neighbor", "WP");
  perfHandle_[static_cast<std::size_t>(ComputePart::DRNeighbor)] =
      run_.flopCounter().addMetric("neighbor-dr", "DR");
}

template <typename Cfg>
void WaveCluster<Cfg>::addReceiverCluster(kernels::ReceiverCluster* receiverCluster) {
  receiverCluster->setNextSampleTime(receiverStartTime_);
  receivers_.push_back(ReceiverSlot{.cluster = receiverCluster, .nextTime = receiverStartTime_});
}

template <typename Cfg>
bool WaveCluster<Cfg>::planPrediction(const StepParams& params) {
  bool outputs = false;
  if (!runTimeOutputs_) {
    for (auto& slot : receivers_) {
      slot.sampling = slot.cluster->planSampling(slot.nextTime, params.time, params.timeStepSize);
      slot.nextTime = slot.sampling.nextTime;
      outputs = outputs || slot.sampling.due;
    }
  }
  incrementPerformanceMetrics(ComputePart::Local);
  return outputs;
}

template <typename Cfg>
void WaveCluster<Cfg>::synchronizeShared(initializer::AllocationPlace place) {
  // the neighbors read the time integrated unknowns or the derivatives in the buffers
  clusterData_->varSynchronizeTo<LTS::Buffers>(place, streamRuntime_.stream());
}

template <typename Cfg>
void WaveCluster<Cfg>::planCorrection(const StepParams& /*params*/) {
  incrementPerformanceMetrics(ComputePart::Neighbor);
  incrementPerformanceMetrics(ComputePart::DRNeighbor);
  // here instead of with the device work, which a replay of a recording skips
  for (auto& stage : cellStages_) {
    stage->countCorrection();
  }
}

template <typename Cfg>
std::function<void()> WaveCluster<Cfg>::hostPrediction(const StepParams& params) {
  // the samples of the receivers as planned for this step
  return [this, params, samplings = plannedSamplings()] {
    writeReceivers(params, samplings);
    computeLocalIntegration(params);
  };
}

template <typename Cfg>
void WaveCluster<Cfg>::devicePrediction(const StepParams& params) {
  writeReceivers(params, plannedSamplings());
  computeLocalIntegrationDevice(params);
}

template <typename Cfg>
void WaveCluster<Cfg>::hostCorrection(const StepParams& params) {
  computeNeighboringIntegration(params);
}

template <typename Cfg>
void WaveCluster<Cfg>::deviceCorrection(const StepParams& params) {
  computeNeighboringIntegrationDevice(params);
}

template <typename Cfg>
void WaveCluster<Cfg>::setTime(double time) {
  CellCluster::setTime(time);
  receiverStartTime_ = time;
  for (auto& slot : receivers_) {
    slot.nextTime = time;
    slot.cluster->setNextSampleTime(time);
  }
}

template <typename Cfg>
void WaveCluster<Cfg>::finalize() {
  cellStages_.clear();
  CellCluster::finalize();
}

template <typename Cfg>
void WaveCluster<Cfg>::finishPhase() {
  CellCluster::finishPhase();
  for (auto& stage : cellStages_) {
    stage->finishPhase();
  }
}

template <typename Cfg>
void WaveCluster<Cfg>::writeReceivers(
    const StepParams& params, const std::vector<kernels::ReceiverCluster::Sampling>& samplings) {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  for (std::size_t i = 0; i < receivers_.size(); ++i) {
    auto* receivers = receivers_[i].cluster;
    if (runTimeOutputs_) {
      receivers->sampleAtRunTime(
          clock_.device(), params.time, params.timeStepSize, executor_, streamRuntime_);
    } else {
      receivers->sample(samplings[i], params.time, params.timeStepSize, executor_, streamRuntime_);
    }
  }
}

template <typename Cfg>
std::vector<kernels::ReceiverCluster::Sampling> WaveCluster<Cfg>::plannedSamplings() const {
  std::vector<kernels::ReceiverCluster::Sampling> samplings;
  samplings.reserve(receivers_.size());
  for (const auto& slot : receivers_) {
    samplings.push_back(slot.sampling);
  }
  return samplings;
}

template <typename Cfg>
void WaveCluster<Cfg>::computeLocalIntegration(const StepParams& params) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  const auto loopBegin = LoopStatistics::begin();

  real* const* stepIntegrals = clusterData_->var<LTS::StepIntegrals>(Cfg());
  real* const* accumulatedIntegrals = clusterData_->var<LTS::AccumulatedIntegrals>(Cfg());
  real* const* derivatives = clusterData_->var<LTS::Derivatives>(Cfg());

  const kernels::LocalTmp<Cfg> tmp(gravitationalAcceleration_);

  const auto timeStepWidth = params.timeStepSize;
  const auto startTime = params.time;
  const auto resetBuffers = params.resetBuffers;
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

  const parallel::runtime::HostRuntime::LoopName loopName(traceLocalIntegration_);
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
void WaveCluster<Cfg>::computeLocalIntegrationDevice(SEISSOL_GPU_PARAM const StepParams& params) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device_.api().putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  const auto loop = deviceLoops_.begin(streamRuntime_.stream());

  auto& dataTable = clusterData_->getConditionalTable();

  kernels::LocalTmp<Cfg> tmp(gravitationalAcceleration_);

  const double timeStepWidth = params.timeStepSize;
  const bool resetBuffers = params.resetBuffers;
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

  const ComputeGraphType graphType =
      resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  const auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, true);
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

        localKernel_.computeBatchedIntegral(dataTable, timeStepWidth, streamRuntime);

        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          const ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
          if (dataTable.contains(key)) {
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
        if (dataTable.contains(key)) {
          auto& entry = dataTable[key];

          if (resetBuffers) {
            device_.algorithms().streamBatchedData(
                const_cast<const real**>(
                    entry.get<real*>(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr()),
                entry.get<real*>(inner_keys::Wp::Id::Buffers)->getDeviceDataPtr(),
                tensor::I<Cfg>::Size,
                entry.get<real*>(inner_keys::Wp::Id::Idofs)->getSize(),
                streamRuntime_.stream());
          } else {
            device_.algorithms().accumulateBatchedData(
                const_cast<const real**>(
                    entry.get<real*>(inner_keys::Wp::Id::Idofs)->getDeviceDataPtr()),
                entry.get<real*>(inner_keys::Wp::Id::Buffers)->getDeviceDataPtr(),
                tensor::I<Cfg>::Size,
                entry.get<real*>(inner_keys::Wp::Id::Idofs)->getSize(),
                streamRuntime_.stream());
          }
        }
      },
      isRecurringTimestep(timeStepWidth));

  // depends on the current time, and therefore cannot be replayed from the graph above; it neither
  // reads nor writes what the graph computes after the local integral
  localKernel_.evaluateBatchedTimeDependentBc(
      dataTable, *clusterData_, clock_.host(), timeStepWidth, streamRuntime_);

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
void WaveCluster<Cfg>::computeNeighboringIntegrationDevice(
    SEISSOL_GPU_PARAM const StepParams& params) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  device_.api().putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  const auto loop = deviceLoops_.begin(streamRuntime_.stream());

  const double timeStepWidth = params.timeStepSize;
  auto& table = clusterData_->getConditionalTable();

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
  const auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime_.runGraph(
      computeGraphKey,
      *clusterData_,
      [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        neighborKernel_.computeBatchedNeighborsIntegral(table, streamRuntime);
      },
      true);

  const bool recurring = isRecurringTimestep(timeStepWidth);
  for (auto& stage : cellStages_) {
    stage->correctDevice(params, recurring, streamRuntime_);
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
void WaveCluster<Cfg>::computeLocalIntegrationFlops() {
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
void WaveCluster<Cfg>::computeNeighborIntegrationFlops() {
  auto& estimateRegular = estimate_[static_cast<int>(ComputePart::Neighbor)];
  auto& estimateDR = estimate_[static_cast<int>(ComputePart::DRNeighbor)];

  estimateRegular = PerformanceEstimate{};
  estimateDR = PerformanceEstimate{};

  const auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  const auto* drMapping = clusterData_->var<LTS::DRMapping>(Cfg());
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    const auto [cellRegular, cellDR] = neighborKernel_.metrics(
        cellInformation[cell].faceTypes, cellInformation[cell].faceRelations, drMapping[cell]);

    estimateRegular += cellRegular;
    estimateDR += cellDR;
  }
}

template <typename Cfg>
void WaveCluster<Cfg>::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();
}

template <typename Cfg>
bool WaveCluster<Cfg>::outputsAhead(long steps) const {
  if (clusterData_->size() == 0 || receivers_.empty() || runTimeOutputs_) {
    return false;
  }
  // the times the host parts of the next steps start with, and the times of the next samples of
  // the outputs
  double time = ct_.correctionTime;
  std::vector<double> nextTimes;
  nextTimes.reserve(receivers_.size());
  for (const auto& slot : receivers_) {
    nextTimes.push_back(slot.nextTime);
  }
  for (long step = 0; step < steps; ++step) {
    const double timeStepSize =
        ct_.timeStepSize(time, ct_.stepsSinceLastSync + step * ct_.timeStepRate, syncTime_);
    for (std::size_t i = 0; i < receivers_.size(); ++i) {
      const auto sampling = receivers_[i].cluster->planSampling(nextTimes[i], time, timeStepSize);
      if (sampling.due) {
        return true;
      }
      nextTimes[i] = sampling.nextTime;
    }
    time += timeStepSize;
  }
  return false;
}

template <typename Cfg>
void WaveCluster<Cfg>::incrementPerformanceMetrics(ComputePart part) {
  run_.flopCounter().incrementMetric(perfHandle_[static_cast<std::size_t>(part)],
                                     estimate_[static_cast<std::size_t>(part)]);
}

template <typename Cfg>
std::vector<parallel::runtime::StreamRuntime*> WaveCluster<Cfg>::runtimes() {
  auto runtimes = CellCluster::runtimes();
  for (const auto& slot : receivers_) {
    if (slot.cluster->sampleRuntime() != nullptr) {
      runtimes.push_back(slot.cluster->sampleRuntime());
    }
  }
  return runtimes;
}

template <typename Cfg>
void WaveCluster<Cfg>::computeNeighboringIntegration(const StepParams& params) {
  const auto clusterSize = clusterData_->size();
  if (clusterSize == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  const auto loopBegin = LoopStatistics::begin();

  const auto* faceNeighbors = clusterData_->var<LTS::FaceNeighbors>();
  const auto* drMapping = clusterData_->var<LTS::DRMapping>(Cfg());
  const auto* cellInformation = clusterData_->var<LTS::CellInformation>();

  // the stages that work on each cell after its correction, and the sums of their counts
  const std::size_t stageCount = cellStages_.size();
  for (auto& stage : cellStages_) {
    stage->beginCorrection(params);
  }
  std::array<std::atomic<std::size_t>, MaxCellStages> stageCounts{};

  const auto timestep = params.timeStepSize;

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs = timeBasis.integrate(
      params.subTimeStart, timestep + params.subTimeStart, params.neighborTimeStepSize);

  const bool configBoundary = !configBoundary_.empty();
  if (configBoundary) {
    configBoundary_.setIntervals(timestep, params.subTimeStart, params.neighborTimeStepSize);
  }

  const parallel::runtime::HostRuntime::LoopName loopName(traceNeighboringIntegration_);
  parallel::runtime::HostRuntime::forRange(clusterSize, [&](std::size_t begin, std::size_t end) {
    std::array<real*, Cell::NumFaces> timeIntegrated{};
    std::array<real*, Cell::NumFaces> faceNeighborsPrefetch{};
    std::array<std::size_t, MaxCellStages> rangeCounts{};

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
                                         data.template get<LTS::Coupling>(),
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

      // the stages work on the cell while its data are still in the cache
      for (std::size_t stage = 0; stage < stageCount; ++stage) {
        rangeCounts[stage] += cellStages_[stage]->correctCell(cell, data);
      }
    }

    for (std::size_t stage = 0; stage < stageCount; ++stage) {
      stageCounts[stage].fetch_add(rangeCounts[stage], std::memory_order_relaxed);
    }
  });

  for (std::size_t stage = 0; stage < stageCount; ++stage) {
    cellStages_[stage]->endCorrection(stageCounts[stage].load(std::memory_order_relaxed));
  }

  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterSize, profilingId_, loopBegin);
}

template <typename Cfg>
void WaveCluster<Cfg>::addCellStage(std::unique_ptr<CellStageBase> stage) {
  if (stage->config() != configIdOf<Cfg>()) {
    logError() << "A cell stage of the configuration" << stage->config()
               << "cannot work on the cells of the configuration" << configIdOf<Cfg>() << ".";
  }
  if (cellStages_.size() == MaxCellStages) {
    logError() << "A cell cluster takes at most" << MaxCellStages << "cell stages.";
  }
  // the configuration is that of the cluster, checked above
  cellStages_.emplace_back(static_cast<CellStage<Cfg>*>(stage.release()));
}

#define SEISSOL_INSTANTIATE(Cfg) template class WaveCluster<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::solver

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
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "Common/Marker.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Precision.h"
#include "Kernels/Receiver.h"
#include "Kernels/Solver.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Parallel/OpenMP.h"
#include "SeisSol.h"
#include "Solver/Settings.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/AbstractAPI.h>
#endif

namespace seissol::time_stepping {

CellCluster::CellCluster(unsigned int clusterId,
                         unsigned int globalClusterId,
                         unsigned int profilingId,
                         const SimulationSettings& settings,
                         HaloType layerType,
                         double maxTimeStepSize,
                         long timeStepRate,
                         bool printProgress,
                         CompoundGlobalData globalData,
                         LTS::Layer* clusterData,
                         seissol::SeisSol& seissolInstance,
                         LoopStatistics* loopStatistics,
                         ActorStateStatistics* actorStateStatistics)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, seissolInstance.executionPlace(clusterData->size())),
      // cluster ids
      settings_(settings), seissolInstance_(seissolInstance), streamRuntime_(4),
      globalDataOnHost_(globalData.onHost), globalDataOnDevice_(globalData.onDevice),
      clusterData_(clusterData),
      sourceCluster_(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      loopStatistics_(loopStatistics), actorStateStatistics_(actorStateStatistics),
      conditionalCounterHost_(1, seissol::memory::Memkind::Standard),
      conditionalCounterDevice_(1,
                                isDeviceOn() ? seissol::memory::Memkind::DeviceGlobalMemory
                                             : seissol::memory::Memkind::Standard),
      layerType_(layerType), printProgress_(printProgress), clusterId_(clusterId),
      globalClusterId_(globalClusterId), profilingId_(profilingId) {
  // assert all pointers are valid
  assert(clusterData_ != nullptr);
  assert(globalDataOnHost_ != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalDataOnDevice_ != nullptr);
  }

  // set timings to zero
  receiverTime_ = 0;

  spacetimeKernel_.setGlobalData(globalData);
  timeKernel_.setGlobalData(globalData);
  localKernel_.setGlobalData(globalData);
  localKernel_.setInitConds(&seissolInstance_.memoryManager().initialConditions());
  localKernel_.setGravitationalAcceleration(seissolInstance_.gravitationSetup().acceleration);
  neighborKernel_.setGlobalData(globalData);

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

void CellCluster::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster_ = std::move(sourceCluster);
}

void CellCluster::writeReceivers(const StepParams& params) {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster_ != nullptr) {
    receiverTime_ = receiverCluster_->calcReceivers(
        receiverTime_, params.time, params.timeStepSize, executor_, streamRuntime_);
  }
}

void CellCluster::computeSources(const StepParams& params) {
#ifdef ACL_DEVICE
  device_.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
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
    loopStatistics_->begin(regionComputePointSources_);
    pointSourceCluster->addTimeIntegratedPointSources(
        params.time, params.time + params.timeStepSize, streamRuntime_);
    loopStatistics_->end(regionComputePointSources_, pointSourceCluster->size(), profilingId_);
  }
#ifdef ACL_DEVICE
  device_.api->popLastProfilingMark();
#endif
}

void CellCluster::computeLocalIntegration(const StepParams& params) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeLocalIntegration_);

  // local integration buffer
  alignas(Alignment) real integrationBuffer[kernels::Solver::IntegralsSize]{};

  // pointer for the call of the ADER-function
  real* bufferPointer = nullptr;

  real* const* stepIntegrals = clusterData_->var<LTS::StepIntegrals>();
  real* const* accumulatedIntegrals = clusterData_->var<LTS::AccumulatedIntegrals>();
  real* const* derivatives = clusterData_->var<LTS::Derivatives>();

  kernels::LocalTmp tmp(seissolInstance_.gravitationSetup().acceleration);

  const auto timeStepWidth = params.timeStepSize;
  const auto startTime = params.time;
  const auto resetBuffers = params.resetBuffers;
  const auto timeBasis = seissol::kernels::timeBasis();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

#pragma omp parallel for private(bufferPointer, integrationBuffer),                                \
    firstprivate(tmp) schedule(static)
  for (std::size_t cell = 0; cell < clusterData_->size(); cell++) {
    auto data = clusterData_->cellRef(cell);

    if (data.get<LTS::CellInformation>().ltsSetup.hasBuffer(BufferType::StepIntegrals)) {
      // assert presence of the buffer
      assert(stepIntegrals[cell] != nullptr);

      bufferPointer = stepIntegrals[cell];
    } else {
      // work on local buffer
      bufferPointer = integrationBuffer;
    }

    spacetimeKernel_.computeAder(
        integrationCoeffs.data(), timeStepWidth, data, tmp, bufferPointer, derivatives[cell], true);

    // Compute local integrals (including local boundary conditions)
    localKernel_.computeIntegral(bufferPointer, data, tmp, startTime, timeStepWidth);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      auto& curFaceDisplacements = data.get<LTS::FaceDisplacements>()[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr &&
          data.get<LTS::CellInformation>().faceTypes[face] != FaceType::FreeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = globalDataOnHost_->v3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
        addVelocityKrnl.faceDisplacement = data.get<LTS::FaceDisplacements>()[face];
        addVelocityKrnl.I = bufferPointer;
        addVelocityKrnl.execute(face);
      }
    }

    // We've used a step integral so far -> accumulate update if needed.
    if (data.get<LTS::CellInformation>().ltsSetup.hasBuffer(BufferType::AccumulatedIntegrals)) {
      assert(accumulatedIntegrals[cell] != nullptr);

      if (resetBuffers) {
        std::memcpy(accumulatedIntegrals[cell],
                    bufferPointer,
                    kernels::Solver::IntegralsSize * sizeof(real));
      } else {
#pragma omp simd
        for (std::size_t dof = 0; dof < kernels::Solver::IntegralsSize; ++dof) {
          accumulatedIntegrals[cell][dof] += bufferPointer[dof];
        }
      }
    }
  }

  loopStatistics_->end(regionComputeLocalIntegration_, clusterData_->size(), profilingId_);
}

void CellCluster::computeLocalIntegrationDevice(SEISSOL_GPU_PARAM const StepParams& params) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device_.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  loopStatistics_->begin(regionComputeLocalIntegration_);

  auto& dataTable = clusterData_->getConditionalTable<inner_keys::Wp>();
  auto& materialTable = clusterData_->getConditionalTable<inner_keys::Material>();
  auto& indicesTable = clusterData_->getConditionalTable<inner_keys::Indices>();

  kernels::LocalTmp tmp(seissolInstance_.gravitationSetup().acceleration);

  const double timeStepWidth = params.timeStepSize;
  const bool resetBuffers = params.resetBuffers;
  const auto timeBasis = seissol::kernels::timeBasis();
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
                                            materialTable,
                                            true,
                                            streamRuntime);

        localKernel_.computeBatchedIntegral(
            dataTable, materialTable, indicesTable, timeStepWidth, streamRuntime);

        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          const ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
          if (dataTable.find(key) != dataTable.end()) {
            auto& entry = dataTable[key];
            // NOTE: integrated velocities have been computed implicitly, i.e
            // it is 6th, 7the and 8th columns of integrated dofs

            kernel::gpu_addVelocity displacementKrnl;
            displacementKrnl.faceDisplacement =
                entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
            displacementKrnl.integratedVelocities = const_cast<const real**>(
                entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
            displacementKrnl.V3mTo2nFace = globalDataOnDevice_->v3mTo2nFace;

            // Note: this kernel doesn't require tmp. memory
            displacementKrnl.numElements =
                entry.get(inner_keys::Wp::Id::FaceDisplacement)->getSize();
            displacementKrnl.streamPtr = streamRuntime_.stream();
            displacementKrnl.execute(face);
          }
        }

        const ConditionalKey key =
            ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
        if (dataTable.find(key) != dataTable.end()) {
          auto& entry = dataTable[key];

          if (resetBuffers) {
            device_.algorithms.streamBatchedData(
                const_cast<const real**>(
                    (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime_.stream());
          } else {
            device_.algorithms.accumulateBatchedData(
                const_cast<const real**>(
                    (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime_.stream());
          }
        }
      });

  // depends on the current time, and therefore cannot be replayed from the graph above; it neither
  // reads nor writes what the graph computes after the local integral
  localKernel_.evaluateBatchedTimeDependentBc(
      dataTable, indicesTable, *clusterData_, params.time, timeStepWidth, streamRuntime_);

  loopStatistics_->end(regionComputeLocalIntegration_, clusterData_->size(), profilingId_);
  device_.api->popLastProfilingMark();
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

void CellCluster::computeNeighboringIntegration(const StepParams& params) {
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

void CellCluster::computeNeighboringIntegrationDevice(SEISSOL_GPU_PARAM const StepParams& params) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  device_.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  loopStatistics_->begin(regionComputeNeighboringIntegration_);

  const double timeStepWidth = params.timeStepSize;
  auto& table = clusterData_->getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
  const auto subtimeCoeffs = timeBasis.integrate(
      params.subTimeStart, timeStepWidth + params.subTimeStart, params.neighborTimeStepSize);

  seissol::kernels::TimeCommon::computeBatchedIntegrals(
      timeKernel_, timeCoeffs.data(), subtimeCoeffs.data(), table, streamRuntime_);

  const ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime_.runGraph(computeGraphKey,
                          *clusterData_,
                          [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
                            neighborKernel_.computeBatchedNeighborsIntegral(table, streamRuntime);
                          });

  if (settings_.plasticity) {
    auto plasticityGraphKey = initializer::GraphKey(ComputeGraphType::Plasticity, timeStepWidth);
    auto* plasticity =
        clusterData_->var<LTS::Plasticity>(seissol::initializer::AllocationPlace::Device);
    auto* isAdjustableVector =
        clusterData_->var<LTS::FlagScratch>(seissol::initializer::AllocationPlace::Device);
    streamRuntime_.runGraph(plasticityGraphKey,
                            *clusterData_,
                            [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
                              seissol::kernels::Plasticity::computePlasticityBatched(
                                  timeStepWidth,
                                  seissolInstance_.parameters().model.tv,
                                  globalDataOnDevice_,
                                  table,
                                  plasticity,
                                  conditionalCounterDevice_.data(),
                                  isAdjustableVector,
                                  streamRuntime);
                            });

    seissolInstance_.flopCounter().incrementMetric(
        perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityCheck)],
        estimate_[static_cast<std::size_t>(ComputePart::PlasticityCheck)] * numPlasticCells_);
  }

  if (settings_.integrate) {
    ConditionalKey key = ConditionalKey(*KernelNames::Time);
    if (table.find(key) != table.end()) {
      auto entry = table.at(key);
      device_.algorithms.accumulateBatchedData(
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr()),
          (entry.get(inner_keys::Wp::Id::Integrals))->getDeviceDataPtr(),
          tensor::Q::Size,
          (entry.get(inner_keys::Wp::Id::Dofs))->getSize(),
          streamRuntime_.stream());
    }
  }

  device_.api->popLastProfilingMark();
  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterData_->size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

void CellCluster::computeLocalIntegrationFlops() {
  auto& estimate = estimate_[static_cast<int>(ComputePart::Local)];
  estimate = PerformanceEstimate{};

  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    estimate += spacetimeKernel_.metrics();
    estimate += localKernel_.metrics(cellInformation[cell].faceTypes);

    // Contribution from displacement/integrated displacement
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::FreeSurfaceGravity) {
        estimate +=
            GravitationalFreeSurfaceBc::metrics(face, cellInformation[cell].faceTypes[face]);
      }
    }
  }
}

void CellCluster::computeNeighborIntegrationFlops() {
  auto& estimateRegular = estimate_[static_cast<int>(ComputePart::Neighbor)];
  auto& estimateDR = estimate_[static_cast<int>(ComputePart::DRNeighbor)];

  estimateRegular = PerformanceEstimate{};
  estimateDR = PerformanceEstimate{};

  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* drMapping = clusterData_->var<LTS::DRMapping>();
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    const auto [cellRegular, cellDR] = neighborKernel_.metrics(
        cellInformation[cell].faceTypes, cellInformation[cell].faceRelations, drMapping[cell]);

    estimateRegular += cellRegular;
    estimateDR += cellDR;
  }
}

void CellCluster::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();

  const auto [check, yield] = seissol::kernels::Plasticity::metrics();
  estimate_[static_cast<int>(ComputePart::PlasticityCheck)] = check;
  estimate_[static_cast<int>(ComputePart::PlasticityYield)] = yield;
}

ActResult CellCluster::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

void CellCluster::handleNeighborPrediction(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
void CellCluster::handleNeighborCorrection(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
void CellCluster::predict() {
  assert(state_ == ActorState::Corrected);
  if (clusterData_->size() == 0) {
    return;
  }

  const auto params = stepParams();

  writeReceivers(params);

  if (executor_ == Executor::Device) {
    computeLocalIntegrationDevice(params);
  } else {
    computeLocalIntegration(params);
  }

  computeSources(params);

  incrementPerformanceMetrics(ComputePart::Local);

  if (hasDifferentExecutorNeighbor()) {
    auto other = executor_ == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                               : seissol::initializer::AllocationPlace::Device;
    clusterData_->varSynchronizeTo<LTS::Buffers>(other, streamRuntime_.stream());
  }

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

void CellCluster::correct() {
  assert(state_ == ActorState::Predicted);
  const auto params = stepParams();

  if (executor_ == Executor::Device) {
    computeNeighboringIntegrationDevice(params);
  } else {
    computeNeighboringIntegration(params);
  }

  incrementPerformanceMetrics(ComputePart::Neighbor);
  incrementPerformanceMetrics(ComputePart::DRNeighbor);

  if (printProgress_) {

    const auto nextCorrectionSteps = ct_.nextCorrectionSteps();
    if (((nextCorrectionSteps / timeStepRate_) % 100) == 0) {
      const auto nextCorrectionTime = ct_.nextCorrectionTime(syncTime_);
      streamRuntime_.enqueueHost([nextCorrectionSteps, nextCorrectionTime]() {
        logInfo() << "Max cluster / LTS cycle updates since sync: " << nextCorrectionSteps
                  << " at time " << nextCorrectionTime;
      });
    }
  }

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

void* CellCluster::recordActionEvent() { return streamRuntime_.eventRecord(); }

void CellCluster::waitForEvent(SEISSOL_GPU_PARAM void* event) {
#ifdef ACL_DEVICE
  if (executor_ == Executor::Host) {
    // the host kernels run right away, so the host has to wait
    device_.api->syncEventWithHost(event);
  } else {
    streamRuntime_.eventSync(event);
  }
#endif
}

void CellCluster::incrementPerformanceMetrics(ComputePart part) {
  seissolInstance_.flopCounter().incrementMetric(perfHandle_[static_cast<std::size_t>(part)],
                                                 estimate_[static_cast<std::size_t>(part)]);
}

unsigned int CellCluster::getClusterId() const { return clusterId_; }

std::size_t CellCluster::layerId() const { return clusterData_->id(); }

unsigned int CellCluster::getGlobalClusterId() const { return globalClusterId_; }

HaloType CellCluster::getLayerType() const { return layerType_; }
void CellCluster::setTime(double time) {
  AbstractTimeCluster::setTime(time);
  this->receiverTime_ = time;
}

void CellCluster::finalize() {
  sourceCluster_.host.reset(nullptr);
  sourceCluster_.device.reset(nullptr);
  streamRuntime_.dispose();

  logDebug() << "#(time steps):" << numberOfTimeSteps_;
}

template <bool UsePlasticity, bool IntegrateOutput>
void CellCluster::computeNeighboringIntegrationImplementation(const StepParams& params) {
  const auto clusterSize = clusterData_->size();
  if (clusterSize == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeNeighboringIntegration_);

  const auto* faceNeighbors = clusterData_->var<LTS::FaceNeighbors>();
  const auto* drMapping = clusterData_->var<LTS::DRMapping>();
  const auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* plasticity = clusterData_->var<LTS::Plasticity>();
  auto* pstrain = clusterData_->var<LTS::PStrain>();

  // NOLINTNEXTLINE
  std::size_t numberOfTetsWithPlasticYielding = 0;

  std::array<real*, Cell::NumFaces> timeIntegrated{};
  std::array<real*, Cell::NumFaces> faceNeighborsPrefetch{};

  const auto tV = seissolInstance_.parameters().model.tv;

  const auto timestep = params.timeStepSize;
  const auto oneMinusIntegratingFactor =
      seissol::kernels::Plasticity::computeRelaxTime(tV, timestep);

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs = timeBasis.integrate(
      params.subTimeStart, timestep + params.subTimeStart, params.neighborTimeStepSize);

#pragma omp parallel for schedule(static) default(none) private(timeIntegrated,                    \
                                                                    faceNeighborsPrefetch)         \
    shared(oneMinusIntegratingFactor,                                                              \
               cellInformation,                                                                    \
               faceNeighbors,                                                                      \
               pstrain,                                                                            \
               plasticity,                                                                         \
               drMapping,                                                                          \
               tV,                                                                                 \
               timeCoeffs,                                                                         \
               subtimeCoeffs,                                                                      \
               clusterData_,                                                                       \
               timestep,                                                                           \
               clusterSize) reduction(+ : numberOfTetsWithPlasticYielding)
  for (std::size_t cell = 0; cell < clusterSize; cell++) {
    auto data = clusterData_->cellRef(cell);

    std::array<real*, Cell::NumFaces> integrationBuffers{};
    for (std::size_t i = 0; i < Cell::NumFaces; ++i) {
      integrationBuffers[i] =
          &globalDataOnHost_->integrationBufferLTS[(OpenMP::threadId() * Cell::NumFaces + i) *
                                                   kernels::Solver::IntegralsSize];
    }

    seissol::kernels::TimeCommon::computeIntegrals(timeKernel_,
                                                   data.get<LTS::CellInformation>().ltsSetup,
                                                   data.get<LTS::CellInformation>().faceTypes,
                                                   timeCoeffs.data(),
                                                   subtimeCoeffs.data(),
                                                   faceNeighbors[cell],
                                                   integrationBuffers,
                                                   timeIntegrated);

    faceNeighborsPrefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                   ? faceNeighbors[cell][1]
                                   : drMapping[cell][1].godunov;
    faceNeighborsPrefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                   ? faceNeighbors[cell][2]
                                   : drMapping[cell][2].godunov;
    faceNeighborsPrefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                   ? faceNeighbors[cell][3]
                                   : drMapping[cell][3].godunov;

    // fourth face's prefetches
    if (cell + 1 < clusterSize) {
      faceNeighborsPrefetch[3] =
          (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
              ? faceNeighbors[cell + 1][0]
              : drMapping[cell + 1][0].godunov;
    } else {
      faceNeighborsPrefetch[3] = faceNeighbors[cell][3];
    }

    neighborKernel_.computeNeighborsIntegral(data, timeIntegrated, faceNeighborsPrefetch);

    if constexpr (UsePlasticity) {
      if (data.get<LTS::CellInformation>().plasticityEnabled) {
        numberOfTetsWithPlasticYielding +=
            seissol::kernels::Plasticity::computePlasticity(oneMinusIntegratingFactor,
                                                            timestep,
                                                            tV,
                                                            globalDataOnHost_,
                                                            &plasticity[cell],
                                                            data.get<LTS::Dofs>(),
                                                            pstrain[cell]);
      }
    }
    if constexpr (IntegrateOutput) {
      auto* __restrict integral = data.get<LTS::Integrals>();
      const auto* __restrict dofs = data.get<LTS::Dofs>();

// only first-order time integration for the output here
#pragma omp simd
      for (std::size_t dof = 0; dof < tensor::Q::size(); ++dof) {
        integral[dof] += timestep * dofs[dof];
      }
    }
  }

  if constexpr (UsePlasticity) {
    conditionalCounterHost_[0] += numberOfTetsWithPlasticYielding;
    seissolInstance_.flopCounter().incrementMetric(
        perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityCheck)],
        estimate_[static_cast<std::size_t>(ComputePart::PlasticityCheck)] * numPlasticCells_);
  }

  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterSize, profilingId_);
}

void CellCluster::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor_ == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor_ == Executor::Host)) {
      clusterData_->synchronizeTo(place, stream);
    }
  }
}

void CellCluster::finishPhase() {
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

std::string CellCluster::description() const {
  const auto identifier = clusterData_->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  return "compute-" + haloStr;
}

} // namespace seissol::time_stepping

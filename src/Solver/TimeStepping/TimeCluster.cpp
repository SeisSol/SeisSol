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

#include "TimeCluster.h"

#include "Alignment.h"
#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "Common/Marker.h"
#include "Config.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/MemoryManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Interface.h"
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Receiver.h"
#include "Kernels/Solver.h"
#include "Kernels/SolverSelector.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Numerical/Quadrature.h"
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

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"

#include <Device/AbstractAPI.h>
#endif

namespace seissol::time_stepping {

template <typename Cfg>
TimeCluster<Cfg>::TimeCluster(
    unsigned int clusterId,
    unsigned int globalClusterId,
    unsigned int profilingId,
    const SimulationSettings& settings,
    HaloType layerType,
    double maxTimeStepSize,
    long timeStepRate,
    bool printProgress,
    DynamicRuptureScheduler* dynamicRuptureScheduler,
    CompoundGlobalData<Cfg> globalData,
    LTS::Layer* clusterData,
    DynamicRupture::Layer* dynRupInteriorData,
    DynamicRupture::Layer* dynRupCopyData,
    const seissol::dr::friction_law::FrictionSolverFactory& frictionSolverFactory,
    const seissol::dr::friction_law::FrictionSolverFactory& frictionSolverFactoryDevice,
    dr::output::OutputManager* faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    ActorStateStatistics* actorStateStatistics)
    : TimeClusterInterface(
          maxTimeStepSize, timeStepRate, seissolInstance.executionPlace(clusterData->size())),
      // cluster ids
      settings_(settings), seissolInstance_(seissolInstance),
      configBoundary_(seissolInstance.parameters().model.configs()), streamRuntime_(4),
      globalData_(globalData), clusterData_(clusterData),
      // global data
      dynRupInteriorData_(dynRupInteriorData), dynRupCopyData_(dynRupCopyData),
      frictionSolver_(dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactory)),
      frictionSolverDevice_(dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactoryDevice)),
      frictionSolverCopy_(dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactory)),
      frictionSolverCopyDevice_(
          dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactoryDevice)),
      faultOutputManager_(faultOutputManager),
      sourceCluster_(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      loopStatistics_(loopStatistics), actorStateStatistics_(actorStateStatistics),
      conditionalCounterHost_(1, seissol::memory::Memkind::Standard),
      conditionalCounterDevice_(1,
                                isDeviceOn() ? seissol::memory::Memkind::DeviceGlobalMemory
                                             : seissol::memory::Memkind::Standard),
      layerType_(layerType), printProgress_(printProgress), clusterId_(clusterId),
      globalClusterId_(globalClusterId), profilingId_(profilingId),
      dynamicRuptureScheduler_(dynamicRuptureScheduler) {
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
  dynamicRuptureKernel_.setGlobalData(globalData);
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

  frictionSolver_->allocateAuxiliaryMemory(globalData_.onHost);
  frictionSolverCopy_->allocateAuxiliaryMemory(globalData_.onHost);
  if constexpr (seissol::isDeviceOn()) {
    frictionSolverDevice_->allocateAuxiliaryMemory(globalData_.onDevice);
    frictionSolverCopyDevice_->allocateAuxiliaryMemory(globalData_.onDevice);
  }

  frictionSolver_->setupLayer(*dynRupInteriorData, streamRuntime_);
  frictionSolverCopy_->setupLayer(*dynRupCopyData, streamRuntime_);
  if constexpr (seissol::isDeviceOn()) {
    frictionSolverDevice_->setupLayer(*dynRupInteriorData, streamRuntime_);
    frictionSolverCopyDevice_->setupLayer(*dynRupCopyData, streamRuntime_);
  }
  streamRuntime_.wait();

  computeFlops();

  regionComputeLocalIntegration_ = loopStatistics_->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration_ =
      loopStatistics_->getRegion("computeNeighboringIntegration");
  regionComputeDynamicRupture_ = loopStatistics_->getRegion("computeDynamicRupture");
  regionComputePointSources_ = loopStatistics_->getRegion("computePointSources");

  conditionalCounterHost_[0] = 0;
  conditionalCounterDevice_.copyFrom(conditionalCounterHost_);

  perfHandle_[static_cast<std::size_t>(ComputePart::Local)] =
      seissolInstance.flopCounter().addMetric("local", "WP");
  perfHandle_[static_cast<std::size_t>(ComputePart::Neighbor)] =
      seissolInstance.flopCounter().addMetric("neighbor", "WP");
  perfHandle_[static_cast<std::size_t>(ComputePart::DRNeighbor)] =
      seissolInstance.flopCounter().addMetric("neighbor-dr", "DR");
  perfHandle_[static_cast<std::size_t>(ComputePart::DRFrictionLawInterior)] =
      seissolInstance.flopCounter().addMetric("dr-frictionlaw-interior", "DR");
  perfHandle_[static_cast<std::size_t>(ComputePart::DRFrictionLawCopy)] =
      seissolInstance.flopCounter().addMetric("dr-frictionlaw-copy", "DR");
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
void TimeCluster<Cfg>::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster_ = std::move(sourceCluster);
}

template <typename Cfg>
void TimeCluster<Cfg>::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster_ != nullptr) {
    receiverTime_ = receiverCluster_->calcReceivers(
        receiverTime_, ct_.correctionTime, timeStepSize(), executor_, streamRuntime_);
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::computeSources() {
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
    loopStatistics_->begin(regionComputePointSources_);
    const auto timeStepSizeLocal = timeStepSize();
    pointSourceCluster->addTimeIntegratedPointSources(
        ct_.correctionTime, ct_.correctionTime + timeStepSizeLocal, streamRuntime_);
    loopStatistics_->end(regionComputePointSources_, pointSourceCluster->size(), profilingId_);
  }
#ifdef ACL_DEVICE
  device_.api().popLastProfilingMark();
#endif
}

template <typename Cfg>
void TimeCluster<Cfg>::computeDynamicRupture(DynamicRupture::Layer& layerData) {
  if (layerData.size() == 0) {
    return;
  }
  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(
      myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON)

  loopStatistics_->begin(regionComputeDynamicRupture_);

  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();
  const DRGodunovData<Cfg>* godunovData = layerData.var<DynamicRupture::GodunovData>(Cfg());
  real* const* timeDerivativePlus = layerData.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real* const* timeDerivativeMinus = layerData.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  auto* qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(Cfg());
  auto* qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(Cfg());

  const auto timestep = timeStepSize();

  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

  const auto pointsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, timestep);
  const auto frictionTime =
      seissol::dr::friction_law::FrictionSolver::computeDeltaT<Cfg>(timePoints);

#pragma omp parallel
  {
    LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation");
  }

#pragma omp parallel for schedule(static)
  for (std::size_t face = 0; face < layerData.size(); ++face) {
    const std::size_t prefetchFace = (face + 1 < layerData.size()) ? face + 1 : face;
    dynamicRuptureKernel_.spaceTimeInterpolation(faceInformation[face],
                                                 &godunovData[face],
                                                 timeDerivativePlus[face],
                                                 timeDerivativeMinus[face],
                                                 qInterpolatedPlus[face],
                                                 qInterpolatedMinus[face],
                                                 timeDerivativePlus[prefetchFace],
                                                 timeDerivativeMinus[prefetchFace],
                                                 pointsCollocate.data());
  }
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("computeDynamicRuptureSpaceTimeInterpolation");
    LIKWID_MARKER_START("computeDynamicRuptureFrictionLaw");
  }

  SCOREP_USER_REGION_BEGIN(
      myRegionHandle, "computeDynamicRuptureFrictionLaw", SCOREP_USER_REGION_TYPE_COMMON)
  auto& solver = &layerData == dynRupInteriorData_ ? frictionSolver_ : frictionSolverCopy_;
  solver->evaluate(ct_.correctionTime, frictionTime, timeWeights.data(), streamRuntime_);
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  loopStatistics_->end(regionComputeDynamicRupture_, layerData.size(), profilingId_);
}

template <typename Cfg>
void TimeCluster<Cfg>::computeDynamicRuptureDevice(
    SEISSOL_GPU_PARAM DynamicRupture::Layer& layerData) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeDynamicRupture_);

  if (layerData.size() > 0) {
    // compute space time interpolation part

    const auto timestep = timeStepSize();

    const ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device_.api().putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, timestep);
    auto& table = layerData.getConditionalTable<inner_keys::Dr>();

    const auto [timePoints, timeWeights] =
        seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

    const auto pointsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, timestep);
    const auto frictionTime =
        seissol::dr::friction_law::FrictionSolver::computeDeltaT<Cfg>(timePoints);

    streamRuntime_.runGraph(
        computeGraphKey,
        layerData,
        [&](seissol::parallel::runtime::StreamRuntime& /*streamRuntime*/) {
          dynamicRuptureKernel_.batchedSpaceTimeInterpolation(
              table, pointsCollocate.data(), streamRuntime_);
        },
        isRecurringTimestep(timestep));
    device_.api().popLastProfilingMark();

    auto& solver =
        &layerData == dynRupInteriorData_ ? frictionSolverDevice_ : frictionSolverCopyDevice_;

    device_.api().putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    if (solver->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo<DynamicRupture::QInterpolatedPlus>(
          initializer::AllocationPlace::Host, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::QInterpolatedMinus>(
          initializer::AllocationPlace::Host, streamRuntime_.stream());
      streamRuntime_.wait();
      solver->evaluate(ct_.correctionTime, frictionTime, timeWeights.data(), streamRuntime_);
      layerData.varSynchronizeTo<DynamicRupture::FluxSolverMinus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::FluxSolverPlus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::ImposedStateMinus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::ImposedStatePlus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
    } else {
      solver->evaluate(ct_.correctionTime, frictionTime, timeWeights.data(), streamRuntime_);
    }

    device_.api().popLastProfilingMark();
  }
  loopStatistics_->end(regionComputeDynamicRupture_, layerData.size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif
}

template <typename Cfg>
PerformanceEstimate TimeCluster<Cfg>::computeDynamicRuptureFlops(DynamicRupture::Layer& layerData) {
  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();

  PerformanceEstimate estimate{};

  for (std::size_t face = 0; face < layerData.size(); ++face) {
    estimate += dynamicRuptureKernel_.metrics(faceInformation[face]);
  }

  return estimate;
}

template <typename Cfg>
void TimeCluster<Cfg>::computeLocalIntegration(bool resetBuffers) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeLocalIntegration_);

  // local integration buffer
  alignas(Alignment) real integrationBuffer[kernels::SolverOf<Cfg>::IntegralsSize]{};

  // pointer for the call of the ADER-function
  real* bufferPointer = nullptr;

  real* const* stepIntegrals = clusterData_->var<LTS::StepIntegrals>(Cfg());
  real* const* accumulatedIntegrals = clusterData_->var<LTS::AccumulatedIntegrals>(Cfg());
  real* const* derivatives = clusterData_->var<LTS::Derivatives>(Cfg());

  kernels::LocalTmp<Cfg> tmp(seissolInstance_.gravitationSetup().acceleration);

  const auto timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

#pragma omp parallel for private(bufferPointer, integrationBuffer),                                \
    firstprivate(tmp) schedule(static)
  for (std::size_t cell = 0; cell < clusterData_->size(); cell++) {
    auto data = clusterData_->cellRef<Cfg>(cell);

    if (data.template get<LTS::CellInformation>().ltsSetup.hasBuffer(BufferType::StepIntegrals)) {
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
    localKernel_.computeIntegral(bufferPointer, data, tmp, ct_.correctionTime, timeStepWidth);

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

  loopStatistics_->end(regionComputeLocalIntegration_, clusterData_->size(), profilingId_);
}

template <typename Cfg>
void TimeCluster<Cfg>::computeLocalIntegrationDevice(SEISSOL_GPU_PARAM bool resetBuffers) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device_.api().putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  loopStatistics_->begin(regionComputeLocalIntegration_);

  auto& dataTable = clusterData_->getConditionalTable<inner_keys::Wp>();
  auto& indicesTable = clusterData_->getConditionalTable<inner_keys::Indices>();

  kernels::LocalTmp<Cfg> tmp(seissolInstance_.gravitationSetup().acceleration);

  const double timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

  // The analytical boundary conditions are evaluated in a host function that is handed the
  // current time. A graph keeps the host functions it recorded as they were, so replaying it
  // would evaluate them at the time of the step that recorded it.
  bool timeDependentBc = false;
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    const ConditionalKey analyticalKey(
        *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
    timeDependentBc = timeDependentBc || indicesTable.find(analyticalKey) != indicesTable.end();
  }

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

        localKernel_.evaluateBatchedTimeDependentBc(dataTable,
                                                    indicesTable,
                                                    *clusterData_,
                                                    ct_.correctionTime,
                                                    timeStepWidth,
                                                    streamRuntime);

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
      isRecurringTimestep(timeStepWidth) && !timeDependentBc);

  loopStatistics_->end(regionComputeLocalIntegration_, clusterData_->size(), profilingId_);
  device_.api().popLastProfilingMark();
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

template <typename Cfg>
void TimeCluster<Cfg>::computeNeighboringIntegration(double subTimeStart) {
  if (settings_.integrate) {
    if (settings_.plasticity) {
      computeNeighboringIntegrationImplementation<true, true>(subTimeStart);
    } else {
      computeNeighboringIntegrationImplementation<false, true>(subTimeStart);
    }
  } else {
    if (settings_.plasticity) {
      computeNeighboringIntegrationImplementation<true, false>(subTimeStart);
    } else {
      computeNeighboringIntegrationImplementation<false, false>(subTimeStart);
    }
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::computeNeighboringIntegrationDevice(SEISSOL_GPU_PARAM double subTimeStart) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  device_.api().putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  loopStatistics_->begin(regionComputeNeighboringIntegration_);

  const double timeStepWidth = timeStepSize();
  auto& table = clusterData_->getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timeStepWidth + subTimeStart, neighborTimestep_);

  seissol::kernels::TimeCommon<Cfg>::computeBatchedIntegrals(
      timeKernel_, timeCoeffs.data(), subtimeCoeffs.data(), table, streamRuntime_);
  if (!configBoundary_.empty()) {
    configBoundary_.setIntervals(timeStepWidth, subTimeStart, neighborTimestep_);
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
    streamRuntime_.runGraph(
        plasticityGraphKey,
        *clusterData_,
        [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
          seissol::kernels::Plasticity<Cfg>::computePlasticityBatched(
              timeStepWidth,
              seissolInstance_.parameters().model.tv,
              globalData_.onDevice,
              table,
              plasticity,
              conditionalCounterDevice_.data(),
              isAdjustableVector,
              streamRuntime);
        },
        isRecurringTimestep(timeStepWidth));

    seissolInstance_.flopCounter().incrementMetric(
        perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityCheck)],
        estimate_[static_cast<std::size_t>(ComputePart::PlasticityCheck)] * numPlasticCells_);
  }

  if (settings_.integrate) {
    ConditionalKey key = ConditionalKey(*KernelNames::Time);
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
  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterData_->size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

template <typename Cfg>
void TimeCluster<Cfg>::computeLocalIntegrationFlops() {
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
void TimeCluster<Cfg>::computeNeighborIntegrationFlops() {
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
void TimeCluster<Cfg>::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();
  estimate_[static_cast<int>(ComputePart::DRFrictionLawInterior)] =
      computeDynamicRuptureFlops(*dynRupInteriorData_);
  estimate_[static_cast<int>(ComputePart::DRFrictionLawCopy)] =
      computeDynamicRuptureFlops(*dynRupCopyData_);

  const auto [check, yield] = seissol::kernels::Plasticity<Cfg>::metrics();
  estimate_[static_cast<int>(ComputePart::PlasticityCheck)] = check;
  estimate_[static_cast<int>(ComputePart::PlasticityYield)] = yield;
}

template <typename Cfg>
ActResult TimeCluster<Cfg>::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

template <typename Cfg>
void TimeCluster<Cfg>::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct_.maxTimeStepSize) {
    lastSubTime_ = neighborCluster.ct.correctionTime;
  }
}
template <typename Cfg>
void TimeCluster<Cfg>::handleAdvancedCorrectionTimeMessage(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
template <typename Cfg>
void TimeCluster<Cfg>::predict() {
  assert(state_ == ActorState::Corrected);
  if (clusterData_->size() == 0) {
    return;
  }

  bool resetBuffers = true;
  for (auto& neighbor : neighbors_) {
    if (neighbor.ct.timeStepRate > ct_.timeStepRate &&
        ct_.stepsSinceLastSync > neighbor.ct.stepsSinceLastSync) {
      resetBuffers = false;
    }
  }
  if (ct_.stepsSinceLastSync == 0) {
    resetBuffers = true;
  }

  writeReceivers();

  if (executor_ == Executor::Device) {
    computeLocalIntegrationDevice(resetBuffers);
  } else {
    computeLocalIntegration(resetBuffers);
  }

  computeSources();

  incrementPerformanceMetrics(ComputePart::Local);

  if (hasDifferentExecutorNeighbor()) {
    auto other = executor_ == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                               : seissol::initializer::AllocationPlace::Device;
    clusterData_->varSynchronizeTo<LTS::Buffers>(other, streamRuntime_.stream());
  }

  streamRuntime_.wait();
}

template <typename Cfg>
void TimeCluster<Cfg>::handleDynamicRupture(DynamicRupture::Layer& layerData) {
  if (layerData.size() == 0) {
    return;
  }

  if (executor_ == Executor::Device) {
    computeDynamicRuptureDevice(layerData);
  } else {
    computeDynamicRupture(layerData);
  }

  double time = ct_.correctionTime;

  // repeat the current solution for some times---to match the existing output scheme.
  // maybe replace with just writePickpointOutput(layerId(), time + dt, dt); some day?

  const double meshDt = ct_.getTimeStepSize();
  // the friction law has just evaluated this step up to its end, and that is the state written out
  const double stateTime = ct_.correctionTime + timeStepSize();

  do {
    const auto oldTime = time;
    time += dynamicRuptureScheduler_->getOutputTimestep();
    const auto trueTime = std::min(time, syncTime_);
    const auto trueDt = trueTime - oldTime;
    faultOutputManager_->writePickpointOutput(
        layerData.id(), stateTime, trueTime, trueDt, meshDt, 0, streamRuntime_);

    // write until we've completed the current copy interval, or if we've hit a sync point
  } while (time * (1 + 1e-8) < ct_.correctionTime + ct_.maxTimeStepSize && time < syncTime_);

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor_ == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                               : seissol::initializer::AllocationPlace::Device;
    layerData.varSynchronizeTo<DynamicRupture::FluxSolverMinus>(other, streamRuntime_.stream());
    layerData.varSynchronizeTo<DynamicRupture::FluxSolverPlus>(other, streamRuntime_.stream());
    layerData.varSynchronizeTo<DynamicRupture::ImposedStateMinus>(other, streamRuntime_.stream());
    layerData.varSynchronizeTo<DynamicRupture::ImposedStatePlus>(other, streamRuntime_.stream());
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::correct() {
  assert(state_ == ActorState::Predicted);
  /* Sub start time of width respect to the next cluster; use 0 if not relevant, for example in GTS.
   * LTS requires to evaluate a partial time integration of the derivatives. The point zero in time
   * refers to the derivation of the surrounding time derivatives, which coincides with the last
   * completed time step of the next cluster. The start/end of the time step is the start/end of
   * this clusters time step relative to the zero point.
   *   Example:
   *                                              5 dt
   *   |-----------------------------------------------------------------------------------------|
   * <<< Time stepping of the next cluster (Cn) (5x larger than the current). |                 | |
   * |                 |                 |
   *   |*****************|*****************|+++++++++++++++++|                 |                 |
   * <<< Status of the current cluster. |                 |                 |                 | | |
   *   |-----------------|-----------------|-----------------|-----------------|-----------------|
   * <<< Time stepping of the current cluster (Cc). 0                 dt               2dt 3dt 4dt
   * 5dt
   *
   *   In the example above two clusters are illustrated: Cc and Cn. Cc is the current cluster under
   * consideration and Cn the next cluster with respect to LTS terminology. Cn is currently at time
   * 0 and provided Cc with derivatives valid until 5dt. Cc updated already twice and did its last
   * full update to reach 2dt (== subTimeStart). Next computeNeighboringCopy is called to accomplish
   * the next full update to reach 3dt (+++). Besides working on the buffers of own buffers and
   * those of previous clusters, Cc needs to evaluate the time prediction of Cn in the interval
   * [2dt, 3dt].
   */
  const double subTimeStart = ct_.correctionTime - lastSubTime_;

  // Note, if this is a copy layer actor, we need the FL_Copy and the FL_Int.
  // Otherwise, this is an interior layer actor, and we need only the FL_Int.
  // We need to avoid computing it twice.
  if (dynamicRuptureScheduler_->mayComputeInterior(ct_.stepsSinceStart)) {
    handleDynamicRupture(*dynRupInteriorData_);

    incrementPerformanceMetrics(ComputePart::DRFrictionLawInterior);

    dynamicRuptureScheduler_->setLastCorrectionStepsInterior(ct_.stepsSinceStart);
  }
  if (layerType_ == HaloType::Copy) {
    handleDynamicRupture(*dynRupCopyData_);

    incrementPerformanceMetrics(ComputePart::DRFrictionLawCopy);

    dynamicRuptureScheduler_->setLastCorrectionStepsCopy((ct_.stepsSinceStart));
  }

  if (executor_ == Executor::Device) {
    computeNeighboringIntegrationDevice(subTimeStart);
  } else {
    computeNeighboringIntegration(subTimeStart);
  }

  incrementPerformanceMetrics(ComputePart::Neighbor);
  incrementPerformanceMetrics(ComputePart::DRNeighbor);

  if (printProgress_) {

    const auto nextCorrectionSteps = ct_.nextCorrectionSteps();
    if (((nextCorrectionSteps / timeStepRate_) % 100) == 0) {
      streamRuntime_.enqueueHost([this, nextCorrectionSteps]() {
        logInfo() << "Max cluster / LTS cycle updates since sync: " << nextCorrectionSteps
                  << " at time " << ct_.nextCorrectionTime(syncTime_);
      });
    }
  }

  streamRuntime_.wait();
}

template <typename Cfg>
void TimeCluster<Cfg>::incrementPerformanceMetrics(ComputePart part) {
  seissolInstance_.flopCounter().incrementMetric(perfHandle_[static_cast<std::size_t>(part)],
                                                 estimate_[static_cast<std::size_t>(part)]);
}

template <typename Cfg>
void TimeCluster<Cfg>::reset() {
  AbstractTimeCluster::reset();
  // note: redundant computation, but it needs to be done somewhere
  neighborTimestep_ = timeStepSize();
  for (auto& neighbor : neighbors_) {
    neighborTimestep_ = std::max(neighbor.ct.getTimeStepSize(), neighborTimestep_);
  }
}

template <typename Cfg>
unsigned int TimeCluster<Cfg>::getClusterId() const {
  return clusterId_;
}

template <typename Cfg>
std::size_t TimeCluster<Cfg>::layerId() const {
  return clusterData_->id();
}

template <typename Cfg>
unsigned int TimeCluster<Cfg>::getGlobalClusterId() const {
  return globalClusterId_;
}

template <typename Cfg>
HaloType TimeCluster<Cfg>::getLayerType() const {
  return layerType_;
}
template <typename Cfg>
void TimeCluster<Cfg>::setTime(double time) {
  AbstractTimeCluster::setTime(time);
  this->receiverTime_ = time;
  this->lastSubTime_ = time;
}

template <typename Cfg>
void TimeCluster<Cfg>::finalize() {
  sourceCluster_.host.reset(nullptr);
  sourceCluster_.device.reset(nullptr);
  streamRuntime_.dispose();

  logDebug() << "#(time steps):" << numberOfTimeSteps_;
}

template <typename Cfg>
template <bool UsePlasticity, bool IntegrateOutput>
void TimeCluster<Cfg>::computeNeighboringIntegrationImplementation(double subTimeStart) {
  const auto clusterSize = clusterData_->size();
  if (clusterSize == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeNeighboringIntegration_);

  const auto* faceNeighbors = clusterData_->var<LTS::FaceNeighbors>();
  const auto* drMapping = clusterData_->var<LTS::DRMapping>(Cfg());
  const auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* plasticity = clusterData_->var<LTS::Plasticity>(Cfg());
  auto* pstrain = clusterData_->var<LTS::PStrain>(Cfg());

  // NOLINTNEXTLINE
  std::size_t numberOfTetsWithPlasticYielding = 0;

  std::array<real*, Cell::NumFaces> timeIntegrated{};
  std::array<real*, Cell::NumFaces> faceNeighborsPrefetch{};

  const auto tV = seissolInstance_.parameters().model.tv;

  const auto timestep = timeStepSize();
  const auto oneMinusIntegratingFactor =
      seissol::kernels::Plasticity<Cfg>::computeRelaxTime(tV, timestep);

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timestep + subTimeStart, neighborTimestep_);

  const bool configBoundary = !configBoundary_.empty();
  if (configBoundary) {
    configBoundary_.setIntervals(timestep, subTimeStart, neighborTimestep_);
  }

#pragma omp parallel for schedule(static) default(none) private(timeIntegrated,                    \
                                                                    faceNeighborsPrefetch)         \
    shared(oneMinusIntegratingFactor,                                                              \
               configBoundary,                                                                     \
               cellInformation,                                                                    \
               faceNeighbors,                                                                      \
               pstrain,                                                                            \
               plasticity,                                                                         \
               drMapping,                                                                          \
               subTimeStart,                                                                       \
               tV,                                                                                 \
               timeCoeffs,                                                                         \
               subtimeCoeffs,                                                                      \
               clusterData_,                                                                       \
               timestep,                                                                           \
               clusterSize) reduction(+ : numberOfTetsWithPlasticYielding)
  for (std::size_t cell = 0; cell < clusterSize; cell++) {
    auto data = clusterData_->cellRef<Cfg>(cell);

    // Scratch for the neighbours whose time integral has to be computed here.
    // Written before it is read, so it needs no initialisation; the frame
    // holds it for the whole loop, one copy per thread.
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

    // fourth face's prefetches
    if (cell + 1 < clusterSize) {
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
        numberOfTetsWithPlasticYielding +=
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
    conditionalCounterHost_[0] += numberOfTetsWithPlasticYielding;
    seissolInstance_.flopCounter().incrementMetric(
        perfHandle_[static_cast<std::size_t>(ComputePart::PlasticityCheck)],
        estimate_[static_cast<std::size_t>(ComputePart::PlasticityCheck)] * numPlasticCells_);
  }

  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterSize, profilingId_);
}

template <typename Cfg>
void TimeCluster<Cfg>::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor_ == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor_ == Executor::Host)) {
      clusterData_->synchronizeTo(place, stream);
      if (layerType_ == HaloType::Interior) {
        dynRupInteriorData_->synchronizeTo(place, stream);
      }
      if (layerType_ == HaloType::Copy) {
        dynRupCopyData_->synchronizeTo(place, stream);
      }
    }
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::finishPhase() {
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
std::string TimeCluster<Cfg>::description() const {
  const auto identifier = clusterData_->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  return "compute-" + haloStr;
}

#define SEISSOL_INSTANTIATE(Cfg) template class TimeCluster<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::time_stepping

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
#include "Common/Constants.h"
#include "Common/Executor.h"
#include "Common/Marker.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Interface.h"
#include "Kernels/LinearCK/GravitationalFreeSurfaceBC.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Precision.h"
#include "Kernels/Receiver.h"
#include "Kernels/Solver.h"
#include "Kernels/TimeCommon.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/MemoryAllocator.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/LoopStatistics.h"
#include "Numerical/Quadrature.h"
#include "Parallel/OpenMP.h"
#include "SeisSol.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
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

TimeCluster::TimeCluster(unsigned int clusterId,
                         unsigned int globalClusterId,
                         unsigned int profilingId,
                         bool usePlasticity,
                         HaloType layerType,
                         double maxTimeStepSize,
                         long timeStepRate,
                         bool printProgress,
                         DynamicRuptureScheduler* dynamicRuptureScheduler,
                         CompoundGlobalData globalData,
                         LTS::Layer* clusterData,
                         DynamicRupture::Layer* dynRupInteriorData,
                         DynamicRupture::Layer* dynRupCopyData,
                         seissol::dr::friction_law::FrictionSolver* frictionSolverTemplate,
                         seissol::dr::friction_law::FrictionSolver* frictionSolverTemplateDevice,
                         dr::output::OutputManager* faultOutputManager,
                         seissol::SeisSol& seissolInstance,
                         LoopStatistics* loopStatistics,
                         ActorStateStatistics* actorStateStatistics)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, seissolInstance.executionPlace(clusterData->size())),
      // cluster ids
      usePlasticity_(usePlasticity), seissolInstance_(seissolInstance), streamRuntime_(4),
      globalDataOnHost_(globalData.onHost), globalDataOnDevice_(globalData.onDevice),
      clusterData_(clusterData),
      // global data
      dynRupInteriorData_(dynRupInteriorData), dynRupCopyData_(dynRupCopyData),
      frictionSolver_(frictionSolverTemplate->clone()),
      frictionSolverDevice_(frictionSolverTemplateDevice->clone()),
      frictionSolverCopy_(frictionSolverTemplate->clone()),
      frictionSolverCopyDevice_(frictionSolverTemplateDevice->clone()),
      faultOutputManager_(faultOutputManager),
      sourceCluster_(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      loopStatistics_(loopStatistics), actorStateStatistics_(actorStateStatistics),
      yieldCells_(1,
                  isDeviceOn() ? seissol::memory::Memkind::PinnedMemory
                               : seissol::memory::Memkind::Standard),
      layerType_(layerType), printProgress_(printProgress), clusterId_(clusterId),
      globalClusterId_(globalClusterId), profilingId_(profilingId),
      dynamicRuptureScheduler_(dynamicRuptureScheduler) {
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
  localKernel_.setInitConds(&seissolInstance_.getMemoryManager().getInitialConditions());
  localKernel_.setGravitationalAcceleration(seissolInstance_.getGravitationSetup().acceleration);
  neighborKernel_.setGlobalData(globalData);
  dynamicRuptureKernel_.setGlobalData(globalData);

  frictionSolver_->allocateAuxiliaryMemory(globalDataOnHost_);
  frictionSolverDevice_->allocateAuxiliaryMemory(globalDataOnDevice_);
  frictionSolverCopy_->allocateAuxiliaryMemory(globalDataOnHost_);
  frictionSolverCopyDevice_->allocateAuxiliaryMemory(globalDataOnDevice_);

  frictionSolver_->setupLayer(*dynRupInteriorData, streamRuntime_);
  frictionSolverDevice_->setupLayer(*dynRupInteriorData, streamRuntime_);
  frictionSolverCopy_->setupLayer(*dynRupCopyData, streamRuntime_);
  frictionSolverCopyDevice_->setupLayer(*dynRupCopyData, streamRuntime_);
  streamRuntime_.wait();

  computeFlops();

  regionComputeLocalIntegration_ = loopStatistics_->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration_ =
      loopStatistics_->getRegion("computeNeighboringIntegration");
  regionComputeDynamicRupture_ = loopStatistics_->getRegion("computeDynamicRupture");
  regionComputePointSources_ = loopStatistics_->getRegion("computePointSources");

  yieldCells_[0] = 0;

  const auto* cellInfo = clusterData_->var<LTS::CellInformation>();
  for (std::size_t i = 0; i < clusterData_->size(); ++i) {
    if (cellInfo[i].plasticityEnabled) {
      ++numPlasticCells_;
    }
  }
}

void TimeCluster::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster_ = std::move(sourceCluster);
}

void TimeCluster::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster_ != nullptr) {
    receiverTime_ = receiverCluster_->calcReceivers(
        receiverTime_, ct_.correctionTime, timeStepSize(), executor_, streamRuntime_);
  }
}

std::vector<NeighborCluster>* TimeCluster::getNeighborClusters() { return &neighbors_; }

void TimeCluster::computeSources() {
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
    const auto timeStepSizeLocal = timeStepSize();
    pointSourceCluster->addTimeIntegratedPointSources(
        ct_.correctionTime, ct_.correctionTime + timeStepSizeLocal, streamRuntime_);
    loopStatistics_->end(regionComputePointSources_, pointSourceCluster->size(), profilingId_);
  }
#ifdef ACL_DEVICE
  device_.api->popLastProfilingMark();
#endif
}

void TimeCluster::computeDynamicRupture(DynamicRupture::Layer& layerData) {
  if (layerData.size() == 0) {
    return;
  }
  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(
      myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON)

  loopStatistics_->begin(regionComputeDynamicRupture_);

  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();
  const DRGodunovData* godunovData = layerData.var<DynamicRupture::GodunovData>();
  real* const* timeDerivativePlus = layerData.var<DynamicRupture::TimeDerivativePlus>();
  real* const* timeDerivativeMinus = layerData.var<DynamicRupture::TimeDerivativeMinus>();
  auto* qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>();
  auto* qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>();

  const auto timestep = timeStepSize();

  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, timestep);

  const auto pointsCollocate = seissol::kernels::timeBasis().collocate(timePoints, timestep);
  const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT(timePoints);

#pragma omp parallel
  {
    LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation");
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
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

void TimeCluster::computeDynamicRuptureDevice(SEISSOL_GPU_PARAM DynamicRupture::Layer& layerData) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeDynamicRupture_);

  if (layerData.size() > 0) {
    // compute space time interpolation part

    const auto timestep = timeStepSize();

    const ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device_.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, timestep);
    auto& table = layerData.getConditionalTable<inner_keys::Dr>();

    const auto [timePoints, timeWeights] =
        seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, timestep);

    const auto pointsCollocate = seissol::kernels::timeBasis().collocate(timePoints, timestep);
    const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT(timePoints);

    streamRuntime_.runGraph(computeGraphKey,
                            layerData,
                            [&](seissol::parallel::runtime::StreamRuntime& /*streamRuntime*/) {
                              dynamicRuptureKernel_.batchedSpaceTimeInterpolation(
                                  table, pointsCollocate.data(), streamRuntime_);
                            });
    device_.api->popLastProfilingMark();

    auto& solver =
        &layerData == dynRupInteriorData_ ? frictionSolverDevice_ : frictionSolverCopyDevice_;

    device_.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
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

    device_.api->popLastProfilingMark();
  }
  loopStatistics_->end(regionComputeDynamicRupture_, layerData.size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif
}

void TimeCluster::computeDynamicRuptureFlops(DynamicRupture::Layer& layerData,
                                             std::uint64_t& nonZeroFlops,
                                             std::uint64_t& hardwareFlops) {
  nonZeroFlops = 0;
  hardwareFlops = 0;

  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();

  for (std::size_t face = 0; face < layerData.size(); ++face) {
    std::uint64_t faceNonZeroFlops = 0;
    std::uint64_t faceHardwareFlops = 0;
    dynamicRuptureKernel_.flopsGodunovState(
        faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

void TimeCluster::computeLocalIntegration(bool resetBuffers) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeLocalIntegration_);

  // local integration buffer
  alignas(Alignment) real integrationBuffer[tensor::I::size()];

  // pointer for the call of the ADER-function
  real* bufferPointer = nullptr;

  real* const* buffers = clusterData_->var<LTS::Buffers>();
  real* const* derivatives = clusterData_->var<LTS::Derivatives>();
  const CellMaterialData* materialData = clusterData_->var<LTS::Material>();

  kernels::LocalTmp tmp(seissolInstance_.getGravitationSetup().acceleration);

  const auto timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

#ifdef _OPENMP
#pragma omp parallel for private(bufferPointer, integrationBuffer),                                \
    firstprivate(tmp) schedule(static)
#endif
  for (std::size_t cell = 0; cell < clusterData_->size(); cell++) {
    auto data = clusterData_->cellRef(cell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided =
        data.get<LTS::CellInformation>().ltsSetup.hasBuffers(); // buffers are provided
    const bool resetMyBuffers =
        buffersProvided && (!data.get<LTS::CellInformation>().ltsSetup.accumulateBuffers() ||
                            resetBuffers); // they should be reset

    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[cell] != nullptr);

      bufferPointer = buffers[cell];
    } else {
      // work on local buffer
      bufferPointer = integrationBuffer;
    }

    spacetimeKernel_.computeAder(
        integrationCoeffs.data(), timeStepWidth, data, tmp, bufferPointer, derivatives[cell], true);

    // Compute local integrals (including some boundary conditions)
    const CellBoundaryMapping(*boundaryMapping)[4] = clusterData_->var<LTS::BoundaryMapping>();
    localKernel_.computeIntegral(bufferPointer,
                                 data,
                                 tmp,
                                 &materialData[cell],
                                 &boundaryMapping[cell],
                                 ct_.correctionTime,
                                 timeStepWidth);

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

    // TODO: Integrate this step into the kernel
    // We've used a temporary buffer -> need to accumulate update in
    // shared buffer.
    if (!resetMyBuffers && buffersProvided) {
      assert(buffers[cell] != nullptr);

      for (std::size_t dof = 0; dof < tensor::I::size(); ++dof) {
        buffers[cell][dof] += integrationBuffer[dof];
      }
    }
  }

  loopStatistics_->end(regionComputeLocalIntegration_, clusterData_->size(), profilingId_);
}

void TimeCluster::computeLocalIntegrationDevice(SEISSOL_GPU_PARAM bool resetBuffers) {

#ifdef ACL_DEVICE
  using namespace seissol::recording;

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device_.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  loopStatistics_->begin(regionComputeLocalIntegration_);

  auto& dataTable = clusterData_->getConditionalTable<inner_keys::Wp>();
  auto& materialTable = clusterData_->getConditionalTable<inner_keys::Material>();
  auto& indicesTable = clusterData_->getConditionalTable<inner_keys::Indices>();

  kernels::LocalTmp tmp(seissolInstance_.getGravitationSetup().acceleration);

  const double timeStepWidth = timeStepSize();
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
                                            tmp,
                                            dataTable,
                                            materialTable,
                                            true,
                                            streamRuntime);

        localKernel_.computeBatchedIntegral(
            dataTable, materialTable, indicesTable, timeStepWidth, streamRuntime);

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

  loopStatistics_->end(regionComputeLocalIntegration_, clusterData_->size(), profilingId_);
  device_.api->popLastProfilingMark();
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

void TimeCluster::computeNeighboringIntegration(double subTimeStart) {
  if (usePlasticity_) {
    computeNeighboringIntegrationImplementation<true>(subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(subTimeStart);
  }
}

void TimeCluster::computeNeighboringIntegrationDevice(SEISSOL_GPU_PARAM double subTimeStart) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;

  device_.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  loopStatistics_->begin(regionComputeNeighboringIntegration_);

  const double timeStepWidth = timeStepSize();
  auto& table = clusterData_->getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timeStepWidth + subTimeStart, neighborTimestep_);

  seissol::kernels::TimeCommon::computeBatchedIntegrals(
      timeKernel_, timeCoeffs.data(), subtimeCoeffs.data(), table, streamRuntime_);

  const ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime_.runGraph(computeGraphKey,
                          *clusterData_,
                          [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
                            neighborKernel_.computeBatchedNeighborsIntegral(table, streamRuntime);
                          });

  if (usePlasticity_) {
    auto plasticityGraphKey = initializer::GraphKey(ComputeGraphType::Plasticity, timeStepWidth);
    auto* plasticity =
        clusterData_->var<LTS::Plasticity>(seissol::initializer::AllocationPlace::Device);
    auto* isAdjustableVector =
        clusterData_->var<LTS::FlagScratch>(seissol::initializer::AllocationPlace::Device);
    streamRuntime_.runGraph(plasticityGraphKey,
                            *clusterData_,
                            [&](seissol::parallel::runtime::StreamRuntime& /*streamRuntime*/) {
                              seissol::kernels::Plasticity::computePlasticityBatched(
                                  timeStepWidth,
                                  seissolInstance_.getSeisSolParameters().model.tv,
                                  globalDataOnDevice_,
                                  table,
                                  plasticity,
                                  yieldCells_.data(),
                                  isAdjustableVector,
                                  streamRuntime_);
                            });

    seissolInstance_.flopCounter().incrementNonZeroFlopsPlasticity(
        numPlasticCells_ * accFlopsNonZero_[static_cast<int>(ComputePart::PlasticityCheck)]);
    seissolInstance_.flopCounter().incrementHardwareFlopsPlasticity(
        numPlasticCells_ * accFlopsHardware_[static_cast<int>(ComputePart::PlasticityCheck)]);
  }

  device_.api->popLastProfilingMark();
  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterData_->size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

void TimeCluster::computeLocalIntegrationFlops() {
  auto& flopsNonZero = accFlopsNonZero_[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = accFlopsHardware_[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    std::uint64_t cellNonZero = 0;
    std::uint64_t cellHardware = 0;
    spacetimeKernel_.flopsAder(cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    localKernel_.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    // Contribution from displacement/integrated displacement
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (cellInformation->faceTypes[face] == FaceType::FreeSurfaceGravity) {
        const auto [nonZeroFlopsDisplacement, hardwareFlopsDisplacement] =
            GravitationalFreeSurfaceBc::getFlopsDisplacementFace(
                face, cellInformation[cell].faceTypes[face]);
        flopsNonZero += nonZeroFlopsDisplacement;
        flopsHardware += hardwareFlopsDisplacement;
      }
    }
  }
}

void TimeCluster::computeNeighborIntegrationFlops() {
  auto& flopsNonZero = accFlopsNonZero_[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = accFlopsHardware_[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = accFlopsNonZero_[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = accFlopsHardware_[static_cast<int>(ComputePart::DRNeighbor)];
  flopsNonZero = 0;
  flopsHardware = 0;
  drFlopsNonZero = 0;
  drFlopsHardware = 0;

  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* drMapping = clusterData_->var<LTS::DRMapping>();
  for (std::size_t cell = 0; cell < clusterData_->size(); ++cell) {
    std::uint64_t cellNonZero = 0;
    std::uint64_t cellHardware = 0;
    std::uint64_t cellDRNonZero = 0;
    std::uint64_t cellDRHardware = 0;
    neighborKernel_.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
                                           cellInformation[cell].faceRelations,
                                           drMapping[cell],
                                           cellNonZero,
                                           cellHardware,
                                           cellDRNonZero,
                                           cellDRHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    drFlopsNonZero += cellDRNonZero;
    drFlopsHardware += cellDRHardware;

    /// \todo add lts time integration
    /// \todo add plasticity
  }
}

void TimeCluster::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();
  computeDynamicRuptureFlops(
      *dynRupInteriorData_,
      accFlopsNonZero_[static_cast<int>(ComputePart::DRFrictionLawInterior)],
      accFlopsHardware_[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
  computeDynamicRuptureFlops(*dynRupCopyData_,
                             accFlopsNonZero_[static_cast<int>(ComputePart::DRFrictionLawCopy)],
                             accFlopsHardware_[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
  seissol::kernels::Plasticity::flopsPlasticity(
      accFlopsNonZero_[static_cast<int>(ComputePart::PlasticityCheck)],
      accFlopsHardware_[static_cast<int>(ComputePart::PlasticityCheck)],
      accFlopsNonZero_[static_cast<int>(ComputePart::PlasticityYield)],
      accFlopsHardware_[static_cast<int>(ComputePart::PlasticityYield)]);
}

ActResult TimeCluster::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

void TimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct_.maxTimeStepSize) {
    lastSubTime_ = neighborCluster.ct.correctionTime;
  }
}
void TimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
void TimeCluster::predict() {
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

  seissolInstance_.flopCounter().incrementNonZeroFlopsLocal(
      accFlopsNonZero_[static_cast<int>(ComputePart::Local)]);
  seissolInstance_.flopCounter().incrementHardwareFlopsLocal(
      accFlopsHardware_[static_cast<int>(ComputePart::Local)]);

  if (hasDifferentExecutorNeighbor()) {
    auto other = executor_ == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                               : seissol::initializer::AllocationPlace::Device;
    clusterData_->varSynchronizeTo<LTS::BuffersDerivatives>(other, streamRuntime_.stream());
  }

  streamRuntime_.wait();
}

void TimeCluster::handleDynamicRupture(DynamicRupture::Layer& layerData) {
  if (executor_ == Executor::Device) {
    computeDynamicRuptureDevice(layerData);
  } else {
    computeDynamicRupture(layerData);
  }

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

void TimeCluster::correct() {
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
    if (dynamicRuptureScheduler_->hasDynamicRuptureFaces()) {
      handleDynamicRupture(*dynRupInteriorData_);
      seissolInstance_.flopCounter().incrementNonZeroFlopsDynamicRupture(
          accFlopsNonZero_[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      seissolInstance_.flopCounter().incrementHardwareFlopsDynamicRupture(
          accFlopsHardware_[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
    }
    dynamicRuptureScheduler_->setLastCorrectionStepsInterior(ct_.stepsSinceStart);
  }
  if (layerType_ == HaloType::Copy) {
    if (dynamicRuptureScheduler_->hasDynamicRuptureFaces()) {
      handleDynamicRupture(*dynRupCopyData_);
      seissolInstance_.flopCounter().incrementNonZeroFlopsDynamicRupture(
          accFlopsNonZero_[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      seissolInstance_.flopCounter().incrementHardwareFlopsDynamicRupture(
          accFlopsHardware_[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
    }
    dynamicRuptureScheduler_->setLastCorrectionStepsCopy((ct_.stepsSinceStart));
  }

  if (executor_ == Executor::Device) {
    computeNeighboringIntegrationDevice(subTimeStart);
  } else {
    computeNeighboringIntegration(subTimeStart);
  }

  seissolInstance_.flopCounter().incrementNonZeroFlopsNeighbor(
      accFlopsNonZero_[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance_.flopCounter().incrementHardwareFlopsNeighbor(
      accFlopsHardware_[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance_.flopCounter().incrementNonZeroFlopsDynamicRupture(
      accFlopsNonZero_[static_cast<int>(ComputePart::DRNeighbor)]);
  seissolInstance_.flopCounter().incrementHardwareFlopsDynamicRupture(
      accFlopsHardware_[static_cast<int>(ComputePart::DRNeighbor)]);

  // First cluster calls fault receiver output
  // Call fault output only if both interior and copy parts of DR were computed
  // TODO: Change from iteration based to time based
  if (dynamicRuptureScheduler_->isFirstClusterWithDynamicRuptureFaces() &&
      dynamicRuptureScheduler_->mayComputeFaultOutput(ct_.stepsSinceStart)) {
    faultOutputManager_->writePickpointOutput(
        ct_.correctionTime + timeStepSize(), timeStepSize(), streamRuntime_);
    dynamicRuptureScheduler_->setLastFaultOutput(ct_.stepsSinceStart);
  }

  streamRuntime_.wait();

  // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
  const auto nextCorrectionSteps = ct_.nextCorrectionSteps();
  if (printProgress_ && (((nextCorrectionSteps / timeStepRate_) % 100) == 0)) {
    logInfo() << "#max-updates since sync: " << nextCorrectionSteps << " @ "
              << ct_.nextCorrectionTime(syncTime_);
  }
}

void TimeCluster::reset() {
  AbstractTimeCluster::reset();
  // note: redundant computation, but it needs to be done somewhere
  neighborTimestep_ = timeStepSize();
  for (auto& neighbor : neighbors_) {
    neighborTimestep_ = std::max(neighbor.ct.getTimeStepSize(), neighborTimestep_);
  }
}

void TimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logWarning(true) << "No update since " << timeSinceLastUpdate.count() << "[s] for global cluster "
                   << globalClusterId_ << " with local cluster id " << clusterId_ << " at state "
                   << actorStateToString(state_) << " predTime = " << ct_.predictionTime
                   << " predictionsSinceSync = " << ct_.predictionsSinceLastSync
                   << " corrTime = " << ct_.correctionTime
                   << " correctionsSinceSync = " << ct_.stepsSinceLastSync
                   << " stepsTillSync = " << ct_.stepsUntilSync << " mayPredict = " << mayPredict()
                   << " mayCorrect = " << mayCorrect() << " maySync = " << maySync();
  for (auto& neighbor : neighbors_) {
    logWarning(true) << "Neighbor with rate = " << neighbor.ct.timeStepRate
                     << "PredTime = " << neighbor.ct.predictionTime
                     << "CorrTime = " << neighbor.ct.correctionTime
                     << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
                     << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }
}

unsigned int TimeCluster::getClusterId() const { return clusterId_; }

std::size_t TimeCluster::layerId() const { return clusterData_->id(); }

unsigned int TimeCluster::getGlobalClusterId() const { return globalClusterId_; }

HaloType TimeCluster::getLayerType() const { return layerType_; }
void TimeCluster::setTime(double time) {
  AbstractTimeCluster::setTime(time);
  this->receiverTime_ = time;
  this->lastSubTime_ = time;
}

void TimeCluster::finalize() {
  sourceCluster_.host.reset(nullptr);
  sourceCluster_.device.reset(nullptr);
  streamRuntime_.dispose();

  logDebug() << "#(time steps):" << numberOfTimeSteps_;
}

template <bool UsePlasticity>
void TimeCluster::computeNeighboringIntegrationImplementation(double subTimeStart) {
  const auto clusterSize = clusterData_->size();
  if (clusterSize == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeNeighboringIntegration_);

  auto* faceNeighbors = clusterData_->var<LTS::FaceNeighbors>();
  auto* drMapping = clusterData_->var<LTS::DRMapping>();
  auto* cellInformation = clusterData_->var<LTS::CellInformation>();
  auto* plasticity = clusterData_->var<LTS::Plasticity>();
  auto* pstrain = clusterData_->var<LTS::PStrain>();

  // NOLINTNEXTLINE
  std::size_t numberOfTetsWithPlasticYielding = 0;

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];

  const auto tV = seissolInstance_.getSeisSolParameters().model.tv;

  const auto timestep = timeStepSize();
  const auto oneMinusIntegratingFactor =
      seissol::kernels::Plasticity::computeRelaxTime(tV, timestep);

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timestep + subTimeStart, neighborTimestep_);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) private(timeIntegrated,                    \
                                                                    faceNeighborsPrefetch)         \
    shared(oneMinusIntegratingFactor,                                                              \
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
#endif
  for (std::size_t cell = 0; cell < clusterSize; cell++) {
    auto data = clusterData_->cellRef(cell);

    seissol::kernels::TimeCommon::computeIntegrals(
        timeKernel_,
        data.get<LTS::CellInformation>().ltsSetup,
        data.get<LTS::CellInformation>().faceTypes,
        timeCoeffs.data(),
        subtimeCoeffs.data(),
        faceNeighbors[cell],
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(
            &(globalDataOnHost_->integrationBufferLTS[OpenMP::threadId() * 4 *
                                                      static_cast<size_t>(tensor::I::size())])),
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

    neighborKernel_.computeNeighborsIntegral(
        data, drMapping[cell], timeIntegrated, faceNeighborsPrefetch);

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
#ifdef INTEGRATE_QUANTITIES
    seissolInstance_.postProcessor().integrateQuantities(
        timeStepWidth_, *clusterData_, cell, dofs[cell]);
#endif // INTEGRATE_QUANTITIES
  }

  if constexpr (UsePlasticity) {
    yieldCells_[0] += numberOfTetsWithPlasticYielding;
    seissolInstance_.flopCounter().incrementNonZeroFlopsPlasticity(
        numPlasticCells_ * accFlopsNonZero_[static_cast<int>(ComputePart::PlasticityCheck)]);
    seissolInstance_.flopCounter().incrementHardwareFlopsPlasticity(
        numPlasticCells_ * accFlopsHardware_[static_cast<int>(ComputePart::PlasticityCheck)]);
  }

  loopStatistics_->end(regionComputeNeighboringIntegration_, clusterSize, profilingId_);
}

void TimeCluster::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
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

void TimeCluster::finishPhase() {
  const auto cells = yieldCells_[0];
  seissolInstance_.flopCounter().incrementNonZeroFlopsPlasticity(
      cells * accFlopsNonZero_[static_cast<int>(ComputePart::PlasticityYield)]);
  seissolInstance_.flopCounter().incrementHardwareFlopsPlasticity(
      cells * accFlopsHardware_[static_cast<int>(ComputePart::PlasticityYield)]);
  yieldCells_[0] = 0;
}

std::string TimeCluster::description() const {
  const auto identifier = clusterData_->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  return "compute-" + haloStr;
}

} // namespace seissol::time_stepping

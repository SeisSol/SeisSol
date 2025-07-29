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
#include "Kernels/DynamicRupture.h"
#include "Kernels/Receiver.h"
#include "Kernels/TimeCommon.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "SeisSol.h"
#include "generated_code/kernel.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/Output/OutputManager.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Interface.h>
#include <Kernels/LinearCK/GravitationalFreeSurfaceBC.h>
#include <Kernels/Plasticity.h>
#include <Kernels/PointSourceCluster.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/MemoryAllocator.h>
#include <Memory/Tree/Layer.h>
#include <Monitoring/ActorStateStatistics.h>
#include <Monitoring/LoopStatistics.h>
#include <Numerical/Quadrature.h>
#include <Solver/TimeStepping/AbstractTimeCluster.h>
#include <Solver/TimeStepping/ActorState.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <init.h>
#include <tensor.h>
#include <utility>
#include <utils/logger.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace seissol::time_stepping {

TimeCluster::TimeCluster(unsigned int clusterId,
                         unsigned int globalClusterId,
                         unsigned int profilingId,
                         bool usePlasticity,
                         LayerType layerType,
                         double maxTimeStepSize,
                         long timeStepRate,
                         bool printProgress,
                         DynamicRuptureScheduler* dynamicRuptureScheduler,
                         CompoundGlobalData globalData,
                         seissol::initializer::Layer* clusterData,
                         seissol::initializer::Layer* dynRupInteriorData,
                         seissol::initializer::Layer* dynRupCopyData,
                         seissol::initializer::LTS* lts,
                         seissol::initializer::DynamicRupture* dynRup,
                         seissol::dr::friction_law::FrictionSolver* frictionSolverTemplate,
                         seissol::dr::friction_law::FrictionSolver* frictionSolverTemplateDevice,
                         dr::output::OutputManager* faultOutputManager,
                         seissol::SeisSol& seissolInstance,
                         LoopStatistics* loopStatistics,
                         ActorStateStatistics* actorStateStatistics)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, seissolInstance.executionPlace(clusterData->size())),
      // cluster ids
      usePlasticity(usePlasticity), seissolInstance(seissolInstance),
      globalDataOnHost(globalData.onHost), globalDataOnDevice(globalData.onDevice),
      clusterData(clusterData),
      // global data
      dynRupInteriorData(dynRupInteriorData), dynRupCopyData(dynRupCopyData), lts(lts),
      dynRup(dynRup), frictionSolver(frictionSolverTemplate->clone()),
      frictionSolverDevice(frictionSolverTemplateDevice->clone()),
      frictionSolverCopy(frictionSolverTemplate->clone()),
      frictionSolverCopyDevice(frictionSolverTemplateDevice->clone()),
      faultOutputManager(faultOutputManager),
      sourceCluster(seissol::kernels::PointSourceClusterPair{nullptr, nullptr}),
      // cells
      loopStatistics(loopStatistics), actorStateStatistics(actorStateStatistics),
      yieldCells(1, isDeviceOn() ? seissol::memory::PinnedMemory : seissol::memory::Standard),
      layerType(layerType), printProgress(printProgress), clusterId(clusterId),
      globalClusterId(globalClusterId), profilingId(profilingId),
      dynamicRuptureScheduler(dynamicRuptureScheduler) {
  // assert all pointers are valid
  assert(clusterData != nullptr);
  assert(globalDataOnHost != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalDataOnDevice != nullptr);
  }

  // set timings to zero
  receiverTime = 0;

  spacetimeKernel.setGlobalData(globalData);
  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  localKernel.setInitConds(&seissolInstance.getMemoryManager().getInitialConditions());
  localKernel.setGravitationalAcceleration(seissolInstance.getGravitationSetup().acceleration);
  neighborKernel.setGlobalData(globalData);
  dynamicRuptureKernel.setGlobalData(globalData);

  frictionSolver->allocateAuxiliaryMemory(globalDataOnHost);
  frictionSolverDevice->allocateAuxiliaryMemory(globalDataOnDevice);
  frictionSolverCopy->allocateAuxiliaryMemory(globalDataOnHost);
  frictionSolverCopyDevice->allocateAuxiliaryMemory(globalDataOnDevice);

  frictionSolver->setupLayer(*dynRupInteriorData, dynRup, streamRuntime);
  frictionSolverDevice->setupLayer(*dynRupInteriorData, dynRup, streamRuntime);
  frictionSolverCopy->setupLayer(*dynRupCopyData, dynRup, streamRuntime);
  frictionSolverCopyDevice->setupLayer(*dynRupCopyData, dynRup, streamRuntime);
  streamRuntime.wait();

  computeFlops();

  regionComputeLocalIntegration = loopStatistics->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration = loopStatistics->getRegion("computeNeighboringIntegration");
  regionComputeDynamicRupture = loopStatistics->getRegion("computeDynamicRupture");
  regionComputePointSources = loopStatistics->getRegion("computePointSources");

  yieldCells[0] = 0;
}

void TimeCluster::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster = std::move(sourceCluster);
}

void TimeCluster::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster != nullptr) {
    receiverTime = receiverCluster->calcReceivers(
        receiverTime, ct.correctionTime, timeStepSize(), executor, nullptr);
  }
}

std::vector<NeighborCluster>* TimeCluster::getNeighborClusters() { return &neighbors; }

void TimeCluster::computeSources() {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
  SCOREP_USER_REGION("computeSources", SCOREP_USER_REGION_TYPE_FUNCTION)

  // Return when point sources not initialized. This might happen if there
  // are no point sources on this rank.
  auto* pointSourceCluster = [&]() -> kernels::PointSourceCluster* {
#ifdef ACL_DEVICE
    if (executor == Executor::Device) {
      return sourceCluster.device.get();
    } else {
      return sourceCluster.host.get();
    }
#else
    return sourceCluster.host.get();
#endif
  }();

  if (pointSourceCluster != nullptr) {
    loopStatistics->begin(regionComputePointSources);
    auto timeStepSizeLocal = timeStepSize();
    pointSourceCluster->addTimeIntegratedPointSources(
        ct.correctionTime, ct.correctionTime + timeStepSizeLocal, streamRuntime);
    loopStatistics->end(regionComputePointSources, pointSourceCluster->size(), profilingId);
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
#endif
}

void TimeCluster::computeDynamicRupture(seissol::initializer::Layer& layerData) {
  if (layerData.size() == 0) {
    return;
  }
  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(
      myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON)

  loopStatistics->begin(regionComputeDynamicRupture);

  DRFaceInformation* faceInformation = layerData.var(dynRup->faceInformation);
  DRGodunovData* godunovData = layerData.var(dynRup->godunovData);
  DREnergyOutput* drEnergyOutput = layerData.var(dynRup->drEnergyOutput);
  real** timeDerivativePlus = layerData.var(dynRup->timeDerivativePlus);
  real** timeDerivativeMinus = layerData.var(dynRup->timeDerivativeMinus);
  auto* qInterpolatedPlus = layerData.var(dynRup->qInterpolatedPlus);
  auto* qInterpolatedMinus = layerData.var(dynRup->qInterpolatedMinus);

  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, timeStepSize());
  const auto pointsCollocate = seissol::kernels::timeBasis().collocate(timePoints, timeStepSize());
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
    dynamicRuptureKernel.spaceTimeInterpolation(faceInformation[face],
                                                globalDataOnHost,
                                                &godunovData[face],
                                                &drEnergyOutput[face],
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
  auto& solver = &layerData == dynRupInteriorData ? frictionSolver : frictionSolverCopy;
  solver->evaluate(ct.correctionTime, frictionTime, timeWeights.data(), streamRuntime);
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  loopStatistics->end(regionComputeDynamicRupture, layerData.size(), profilingId);
}

#ifdef ACL_DEVICE
void TimeCluster::computeDynamicRuptureDevice(seissol::initializer::Layer& layerData) {
  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeDynamicRupture);

  if (layerData.size() > 0) {
    // compute space time interpolation part

    const double stepSizeWidth = timeStepSize();
    ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, stepSizeWidth);
    auto& table = layerData.getConditionalTable<inner_keys::Dr>();

    const auto [timePoints, timeWeights] =
        seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, timeStepSize());
    const auto pointsCollocate = seissol::kernels::timeBasis().collocate(timePoints, stepSizeWidth);
    const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT(timePoints);

    streamRuntime.runGraph(
        computeGraphKey, layerData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
          dynamicRuptureKernel.batchedSpaceTimeInterpolation(
              table, pointsCollocate.data(), streamRuntime);
        });
    device.api->popLastProfilingMark();
    if (frictionSolverDevice->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo(
          dynRup->qInterpolatedPlus, initializer::AllocationPlace::Host, streamRuntime.stream());
      layerData.varSynchronizeTo(
          dynRup->qInterpolatedMinus, initializer::AllocationPlace::Host, streamRuntime.stream());
      streamRuntime.wait();
    }

    device.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    auto& solver =
        &layerData == dynRupInteriorData ? frictionSolverDevice : frictionSolverCopyDevice;
    solver->evaluate(ct.correctionTime, frictionTime, timeWeights.data(), streamRuntime);
    device.api->popLastProfilingMark();
    if (frictionSolverDevice->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo(
          dynRup->fluxSolverMinus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo(
          dynRup->fluxSolverPlus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo(
          dynRup->imposedStateMinus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo(
          dynRup->imposedStatePlus, initializer::AllocationPlace::Device, streamRuntime.stream());
    }
    streamRuntime.wait();
  }
  loopStatistics->end(regionComputeDynamicRupture, layerData.size(), profilingId);
}
#endif

void TimeCluster::computeDynamicRuptureFlops(seissol::initializer::Layer& layerData,
                                             std::uint64_t& nonZeroFlops,
                                             std::uint64_t& hardwareFlops) {
  nonZeroFlops = 0;
  hardwareFlops = 0;

  DRFaceInformation* faceInformation = layerData.var(dynRup->faceInformation);

  for (std::size_t face = 0; face < layerData.size(); ++face) {
    std::uint64_t faceNonZeroFlops = 0;
    std::uint64_t faceHardwareFlops = 0;
    dynamicRuptureKernel.flopsGodunovState(
        faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

void TimeCluster::computeLocalIntegration(bool resetBuffers) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeLocalIntegration);

  // local integration buffer
  alignas(Alignment) real integrationBuffer[tensor::I::size()];

  // pointer for the call of the ADER-function
  real* bufferPointer = nullptr;

  real** buffers = clusterData->var(lts->buffers);
  real** derivatives = clusterData->var(lts->derivatives);
  CellMaterialData* materialData = clusterData->var(lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*lts, *clusterData);
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

  const auto timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

#ifdef _OPENMP
#pragma omp parallel for private(bufferPointer, integrationBuffer),                                \
    firstprivate(tmp) schedule(static)
#endif
  for (std::size_t cell = 0; cell < clusterData->size(); cell++) {
    auto data = loader.entry(cell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided =
        (data.cellInformation().ltsSetup >> 8) % 2 == 1; // buffers are provided
    const bool resetMyBuffers =
        buffersProvided &&
        ((data.cellInformation().ltsSetup >> 10) % 2 == 0 || resetBuffers); // they should be reset

    if (resetMyBuffers) {
      // assert presence of the buffer
      assert(buffers[cell] != nullptr);

      bufferPointer = buffers[cell];
    } else {
      // work on local buffer
      bufferPointer = integrationBuffer;
    }

    spacetimeKernel.computeAder(
        integrationCoeffs.data(), timeStepWidth, data, tmp, bufferPointer, derivatives[cell], true);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping(*boundaryMapping)[4] = clusterData->var(lts->boundaryMapping);
    localKernel.computeIntegral(bufferPointer,
                                data,
                                tmp,
                                &materialData[cell],
                                &boundaryMapping[cell],
                                ct.correctionTime,
                                timeStepWidth);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      auto& curFaceDisplacements = data.faceDisplacements()[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr &&
          data.cellInformation().faceTypes[face] != FaceType::FreeSurfaceGravity) {
        kernel::addVelocity addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = globalDataOnHost->v3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity::Values;
        addVelocityKrnl.faceDisplacement = data.faceDisplacements()[face];
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

  loopStatistics->end(regionComputeLocalIntegration, clusterData->size(), profilingId);
}
#ifdef ACL_DEVICE
void TimeCluster::computeLocalIntegrationDevice(bool resetBuffers) {

  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  loopStatistics->begin(regionComputeLocalIntegration);

  auto& dataTable = clusterData->getConditionalTable<inner_keys::Wp>();
  auto& materialTable = clusterData->getConditionalTable<inner_keys::Material>();
  auto& indicesTable = clusterData->getConditionalTable<inner_keys::Indices>();

  kernels::LocalData::Loader loader;
  loader.load(*lts, *clusterData);
  kernels::LocalTmp tmp(seissolInstance.getGravitationSetup().acceleration);

  const double timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

  ComputeGraphType graphType =
      resetBuffers ? ComputeGraphType::AccumulatedVelocities : ComputeGraphType::StreamedVelocities;
  auto computeGraphKey = initializer::GraphKey(graphType, timeStepWidth, true);
  streamRuntime.runGraph(
      computeGraphKey, *clusterData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        spacetimeKernel.computeBatchedAder(integrationCoeffs.data(),
                                           timeStepWidth,
                                           tmp,
                                           dataTable,
                                           materialTable,
                                           true,
                                           streamRuntime);

        localKernel.computeBatchedIntegral(
            dataTable, materialTable, indicesTable, loader, tmp, timeStepWidth, streamRuntime);

        localKernel.evaluateBatchedTimeDependentBc(dataTable,
                                                   indicesTable,
                                                   loader,
                                                   *clusterData,
                                                   *lts,
                                                   ct.correctionTime,
                                                   timeStepWidth,
                                                   streamRuntime);

        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
          if (dataTable.find(key) != dataTable.end()) {
            auto& entry = dataTable[key];
            // NOTE: integrated velocities have been computed implicitly, i.e
            // it is 6th, 7the and 8th columns of integrated dofs

            kernel::gpu_addVelocity displacementKrnl;
            displacementKrnl.faceDisplacement =
                entry.get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
            displacementKrnl.integratedVelocities = const_cast<const real**>(
                entry.get(inner_keys::Wp::Id::Ivelocities)->getDeviceDataPtr());
            displacementKrnl.V3mTo2nFace = globalDataOnDevice->v3mTo2nFace;

            // Note: this kernel doesn't require tmp. memory
            displacementKrnl.numElements =
                entry.get(inner_keys::Wp::Id::FaceDisplacement)->getSize();
            displacementKrnl.streamPtr = streamRuntime.stream();
            displacementKrnl.execute(face);
          }
        }

        ConditionalKey key = ConditionalKey(*KernelNames::Time, *ComputationKind::WithLtsBuffers);
        if (dataTable.find(key) != dataTable.end()) {
          auto& entry = dataTable[key];

          if (resetBuffers) {
            device.algorithms.streamBatchedData(
                const_cast<const real**>(
                    (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime.stream());
          } else {
            device.algorithms.accumulateBatchedData(
                const_cast<const real**>(
                    (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime.stream());
          }
        }
      });

  streamRuntime.wait();

  loopStatistics->end(regionComputeLocalIntegration, clusterData->size(), profilingId);
  device.api->popLastProfilingMark();
}
#endif // ACL_DEVICE

void TimeCluster::computeNeighboringIntegration(double subTimeStart) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(subTimeStart);
  }
}
#ifdef ACL_DEVICE
void TimeCluster::computeNeighboringIntegrationDevice(double subTimeStart) {
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  loopStatistics->begin(regionComputeNeighboringIntegration);

  const double timeStepWidth = timeStepSize();
  auto& table = clusterData->getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timeStepWidth + subTimeStart, neighborTimestep);

  seissol::kernels::TimeCommon::computeBatchedIntegrals(
      timeKernel, timeCoeffs.data(), subtimeCoeffs.data(), table, streamRuntime);

  ComputeGraphType graphType = ComputeGraphType::NeighborIntegral;
  auto computeGraphKey = initializer::GraphKey(graphType);

  streamRuntime.runGraph(
      computeGraphKey, *clusterData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
        neighborKernel.computeBatchedNeighborsIntegral(table, streamRuntime);
      });

  if (usePlasticity) {
    auto plasticityGraphKey = initializer::GraphKey(ComputeGraphType::Plasticity, timeStepWidth);
    auto* plasticity =
        clusterData->var(lts->plasticity, seissol::initializer::AllocationPlace::Device);
    auto* isAdjustableVector =
        clusterData->var(lts->flagScratch, seissol::initializer::AllocationPlace::Device);
    streamRuntime.runGraph(plasticityGraphKey,
                           *clusterData,
                           [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
                             seissol::kernels::Plasticity::computePlasticityBatched(
                                 timeStepWidth,
                                 seissolInstance.getSeisSolParameters().model.tv,
                                 globalDataOnDevice,
                                 table,
                                 plasticity,
                                 yieldCells.data(),
                                 isAdjustableVector,
                                 streamRuntime);
                           });

    seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
        clusterData->size() * accFlopsNonZero[static_cast<int>(ComputePart::PlasticityCheck)]);
    seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
        clusterData->size() * accFlopsHardware[static_cast<int>(ComputePart::PlasticityCheck)]);
  }

  device.api->popLastProfilingMark();
  streamRuntime.wait();
  loopStatistics->end(regionComputeNeighboringIntegration, clusterData->size(), profilingId);
}
#endif // ACL_DEVICE

void TimeCluster::computeLocalIntegrationFlops() {
  auto& flopsNonZero = accFlopsNonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = accFlopsHardware[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = clusterData->var(lts->cellInformation);
  for (std::size_t cell = 0; cell < clusterData->size(); ++cell) {
    std::uint64_t cellNonZero = 0;
    std::uint64_t cellHardware = 0;
    spacetimeKernel.flopsAder(cellNonZero, cellHardware);
    flopsNonZero += cellNonZero;
    flopsHardware += cellHardware;
    localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
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
  auto& flopsNonZero = accFlopsNonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = accFlopsHardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = accFlopsNonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = accFlopsHardware[static_cast<int>(ComputePart::DRNeighbor)];
  flopsNonZero = 0;
  flopsHardware = 0;
  drFlopsNonZero = 0;
  drFlopsHardware = 0;

  auto* cellInformation = clusterData->var(lts->cellInformation);
  auto* drMapping = clusterData->var(lts->drMapping);
  for (std::size_t cell = 0; cell < clusterData->size(); ++cell) {
    std::uint64_t cellNonZero = 0;
    std::uint64_t cellHardware = 0;
    std::uint64_t cellDRNonZero = 0;
    std::uint64_t cellDRHardware = 0;
    neighborKernel.flopsNeighborsIntegral(cellInformation[cell].faceTypes,
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
      *dynRupInteriorData,
      accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)],
      accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
  computeDynamicRuptureFlops(*dynRupCopyData,
                             accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)],
                             accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
  seissol::kernels::Plasticity::flopsPlasticity(
      accFlopsNonZero[static_cast<int>(ComputePart::PlasticityCheck)],
      accFlopsHardware[static_cast<int>(ComputePart::PlasticityCheck)],
      accFlopsNonZero[static_cast<int>(ComputePart::PlasticityYield)],
      accFlopsHardware[static_cast<int>(ComputePart::PlasticityYield)]);
}

ActResult TimeCluster::act() {
  actorStateStatistics->enter(state);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics->enter(state);
  return result;
}

void TimeCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct.maxTimeStepSize) {
    lastSubTime = neighborCluster.ct.correctionTime;
  }
}
void TimeCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}
void TimeCluster::predict() {
  assert(state == ActorState::Corrected);
  if (clusterData->size() == 0) {
    return;
  }

  bool resetBuffers = true;
  for (auto& neighbor : neighbors) {
    if (neighbor.ct.timeStepRate > ct.timeStepRate &&
        ct.stepsSinceLastSync > neighbor.ct.stepsSinceLastSync) {
      resetBuffers = false;
    }
  }
  if (ct.stepsSinceLastSync == 0) {
    resetBuffers = true;
  }

  writeReceivers();
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeLocalIntegrationDevice(resetBuffers);
  } else {
    computeLocalIntegration(resetBuffers);
  }
#else
  computeLocalIntegration(resetBuffers);
#endif
  computeSources();

  seissolInstance.flopCounter().incrementNonZeroFlopsLocal(
      accFlopsNonZero[static_cast<int>(ComputePart::Local)]);
  seissolInstance.flopCounter().incrementHardwareFlopsLocal(
      accFlopsHardware[static_cast<int>(ComputePart::Local)]);
#ifdef ACL_DEVICE
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                              : seissol::initializer::AllocationPlace::Device;
    clusterData->varSynchronizeTo(lts->buffersDerivatives, other, streamRuntime.stream());
    streamRuntime.wait();
  }
#endif
}

void TimeCluster::handleDynamicRupture(initializer::Layer& layerData) {
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeDynamicRuptureDevice(layerData);
  } else {
    computeDynamicRupture(layerData);
  }

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                              : seissol::initializer::AllocationPlace::Device;
    layerData.varSynchronizeTo(dynRup->fluxSolverMinus, other, streamRuntime.stream());
    layerData.varSynchronizeTo(dynRup->fluxSolverPlus, other, streamRuntime.stream());
    layerData.varSynchronizeTo(dynRup->imposedStateMinus, other, streamRuntime.stream());
    layerData.varSynchronizeTo(dynRup->imposedStatePlus, other, streamRuntime.stream());
    streamRuntime.wait();
  }
#else
  computeDynamicRupture(layerData);
#endif
}

void TimeCluster::correct() {
  assert(state == ActorState::Predicted);
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
  const double subTimeStart = ct.correctionTime - lastSubTime;

  // Note, if this is a copy layer actor, we need the FL_Copy and the FL_Int.
  // Otherwise, this is an interior layer actor, and we need only the FL_Int.
  // We need to avoid computing it twice.
  if (dynamicRuptureScheduler->mayComputeInterior(ct.stepsSinceStart)) {
    if (dynamicRuptureScheduler->hasDynamicRuptureFaces()) {
      handleDynamicRupture(*dynRupInteriorData);
      seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(
          accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
      seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(
          accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
    }
    dynamicRuptureScheduler->setLastCorrectionStepsInterior(ct.stepsSinceStart);
  }
  if (layerType == Copy) {
    if (dynamicRuptureScheduler->hasDynamicRuptureFaces()) {
      handleDynamicRupture(*dynRupCopyData);
      seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(
          accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(
          accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
    }
    dynamicRuptureScheduler->setLastCorrectionStepsCopy((ct.stepsSinceStart));
  }

#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeNeighboringIntegrationDevice(subTimeStart);
  } else {
    computeNeighboringIntegration(subTimeStart);
  }
#else
  computeNeighboringIntegration(subTimeStart);
#endif

  seissolInstance.flopCounter().incrementNonZeroFlopsNeighbor(
      accFlopsNonZero[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance.flopCounter().incrementHardwareFlopsNeighbor(
      accFlopsHardware[static_cast<int>(ComputePart::Neighbor)]);
  seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(
      accFlopsNonZero[static_cast<int>(ComputePart::DRNeighbor)]);
  seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(
      accFlopsHardware[static_cast<int>(ComputePart::DRNeighbor)]);

  // First cluster calls fault receiver output
  // Call fault output only if both interior and copy parts of DR were computed
  // TODO: Change from iteration based to time based
  if (dynamicRuptureScheduler->isFirstClusterWithDynamicRuptureFaces() &&
      dynamicRuptureScheduler->mayComputeFaultOutput(ct.stepsSinceStart)) {
    faultOutputManager->writePickpointOutput(ct.correctionTime + timeStepSize(), timeStepSize());
    dynamicRuptureScheduler->setLastFaultOutput(ct.stepsSinceStart);
  }

  // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
  const auto nextCorrectionSteps = ct.nextCorrectionSteps();
  if (printProgress && (((nextCorrectionSteps / timeStepRate) % 100) == 0)) {
    logInfo() << "#max-updates since sync: " << nextCorrectionSteps << " @ "
              << ct.nextCorrectionTime(syncTime);
  }
}

void TimeCluster::reset() {
  AbstractTimeCluster::reset();
  // note: redundant computation, but it needs to be done somewhere
  neighborTimestep = timeStepSize();
  for (auto& neighbor : neighbors) {
    neighborTimestep = std::max(neighbor.ct.getTimeStepSize(), neighborTimestep);
  }
}

void TimeCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logWarning(true) << "No update since " << timeSinceLastUpdate.count() << "[s] for global cluster "
                   << globalClusterId << " with local cluster id " << clusterId << " at state "
                   << actorStateToString(state) << " predTime = " << ct.predictionTime
                   << " predictionsSinceSync = " << ct.predictionsSinceLastSync
                   << " corrTime = " << ct.correctionTime
                   << " correctionsSinceSync = " << ct.stepsSinceLastSync
                   << " stepsTillSync = " << ct.stepsUntilSync << " mayPredict = " << mayPredict()
                   << " mayCorrect = " << mayCorrect() << " maySync = " << maySync();
  for (auto& neighbor : neighbors) {
    logWarning(true) << "Neighbor with rate = " << neighbor.ct.timeStepRate
                     << "PredTime = " << neighbor.ct.predictionTime
                     << "CorrTime = " << neighbor.ct.correctionTime
                     << "predictionsSinceSync = " << neighbor.ct.predictionsSinceLastSync
                     << "correctionsSinceSync = " << neighbor.ct.stepsSinceLastSync;
  }
}

unsigned int TimeCluster::getClusterId() const { return clusterId; }

unsigned int TimeCluster::getGlobalClusterId() const { return globalClusterId; }

LayerType TimeCluster::getLayerType() const { return layerType; }
void TimeCluster::setTime(double time) {
  AbstractTimeCluster::setTime(time);
  this->receiverTime = time;
  this->lastSubTime = time;
}

void TimeCluster::finalize() {
  sourceCluster.host.reset(nullptr);
  sourceCluster.device.reset(nullptr);
  streamRuntime.dispose();

  logDebug() << "#(time steps):" << numberOfTimeSteps;
}

template <bool UsePlasticity>
void TimeCluster::computeNeighboringIntegrationImplementation(double subTimeStart) {
  if (clusterData->size() == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeNeighboringIntegration);

  auto* faceNeighbors = clusterData->var(lts->faceNeighbors);
  auto* drMapping = clusterData->var(lts->drMapping);
  auto* cellInformation = clusterData->var(lts->cellInformation);
  auto* plasticity = clusterData->var(lts->plasticity);
  auto* pstrain = clusterData->var(lts->pstrain);

  // NOLINTNEXTLINE
  std::size_t numberOTetsWithPlasticYielding = 0;

  kernels::NeighborData::Loader loader;
  loader.load(*lts, *clusterData);

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];

  const auto tV = seissolInstance.getSeisSolParameters().model.tv;

  const auto timestep = timeStepSize();
  const auto oneMinusIntegratingFactor =
      seissol::kernels::Plasticity::computeRelaxTime(tV, timeStepSize());

  const auto timeBasis = seissol::kernels::timeBasis();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timestep + subTimeStart, neighborTimestep);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) private(timeIntegrated,                    \
                                                                    faceNeighborsPrefetch)         \
    shared(oneMinusIntegratingFactor,                                                              \
               cellInformation,                                                                    \
               loader,                                                                             \
               faceNeighbors,                                                                      \
               pstrain,                                                                            \
               plasticity,                                                                         \
               drMapping,                                                                          \
               subTimeStart,                                                                       \
               tV,                                                                                 \
               timeCoeffs,                                                                         \
               subtimeCoeffs) reduction(+ : numberOTetsWithPlasticYielding)
#endif
  for (std::size_t cell = 0; cell < clusterData->size(); cell++) {
    auto data = loader.entry(cell);

#ifdef _OPENMP
    const std::size_t threadId = omp_get_thread_num();
#else
    const std::size_t threadId = 0;
#endif
    seissol::kernels::TimeCommon::computeIntegrals(
        timeKernel,
        data.cellInformation().ltsSetup,
        data.cellInformation().faceTypes,
        timeCoeffs.data(),
        subtimeCoeffs.data(),
        faceNeighbors[cell],
        *reinterpret_cast<real(*)[4][tensor::I::size()]>(
            &(globalDataOnHost->integrationBufferLTS[threadId * 4 * tensor::I::size()])),
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
    if (cell + 1 < clusterData->size()) {
      faceNeighborsPrefetch[3] =
          (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
              ? faceNeighbors[cell + 1][0]
              : drMapping[cell + 1][0].godunov;
    } else {
      faceNeighborsPrefetch[3] = faceNeighbors[cell][3];
    }

    neighborKernel.computeNeighborsIntegral(
        data, drMapping[cell], timeIntegrated, faceNeighborsPrefetch);

    if constexpr (UsePlasticity) {
      numberOTetsWithPlasticYielding +=
          seissol::kernels::Plasticity::computePlasticity(oneMinusIntegratingFactor,
                                                          timeStepSize(),
                                                          tV,
                                                          globalDataOnHost,
                                                          &plasticity[cell],
                                                          data.dofs(),
                                                          pstrain[cell]);
    }
#ifdef INTEGRATE_QUANTITIES
    seissolInstance.postProcessor().integrateQuantities(
        m_timeStepWidth, *clusterData, cell, dofs[cell]);
#endif // INTEGRATE_QUANTITIES
  }

  if constexpr (UsePlasticity) {
    yieldCells[0] += numberOTetsWithPlasticYielding;
    seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
        clusterData->size() * accFlopsNonZero[static_cast<int>(ComputePart::PlasticityCheck)]);
    seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
        clusterData->size() * accFlopsHardware[static_cast<int>(ComputePart::PlasticityCheck)]);
  }

  loopStatistics->end(regionComputeNeighboringIntegration, clusterData->size(), profilingId);
}

void TimeCluster::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
#ifdef ACL_DEVICE
  if ((place == initializer::AllocationPlace::Host && executor == Executor::Device) ||
      (place == initializer::AllocationPlace::Device && executor == Executor::Host)) {
    clusterData->synchronizeTo(place, stream);
    if (layerType == Interior) {
      dynRupInteriorData->synchronizeTo(place, stream);
    }
    if (layerType == Copy) {
      dynRupCopyData->synchronizeTo(place, stream);
    }
  }
#endif
}

void TimeCluster::finishPhase() {
  const auto cells = yieldCells[0];
  seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
      cells * accFlopsNonZero[static_cast<int>(ComputePart::PlasticityYield)]);
  seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
      cells * accFlopsHardware[static_cast<int>(ComputePart::PlasticityYield)]);
  yieldCells[0] = 0;
}

} // namespace seissol::time_stepping

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
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Receiver.h"
#include "Kernels/TimeCommon.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Parallel/OpenMP.h"
#include "SeisSol.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <Common/Executor.h>
#include <Config.h>
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/Output/OutputManager.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Interface.h>
#include <Kernels/LinearCK/GravitationalFreeSurfaceBC.h>
#include <Kernels/Plasticity.h>
#include <Kernels/PointSourceCluster.h>
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
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::time_stepping {

template <typename Cfg>
TimeCluster<Cfg>::TimeCluster(
    unsigned int clusterId,
    unsigned int globalClusterId,
    unsigned int profilingId,
    bool usePlasticity,
    HaloType layerType,
    double maxTimeStepSize,
    long timeStepRate,
    bool printProgress,
    DynamicRuptureScheduler* dynamicRuptureScheduler,
    const GlobalData& globalData,
    LTS::Layer* clusterData,
    DynamicRupture::Layer* dynRupInteriorData,
    DynamicRupture::Layer* dynRupCopyData,
    seissol::dr::friction_law::FrictionSolverFactory& frictionSolverTemplate,
    seissol::dr::friction_law::FrictionSolverFactory& frictionSolverTemplateDevice,
    dr::output::OutputManager* faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    ActorStateStatistics* actorStateStatistics)
    : TimeClusterInterface(
          maxTimeStepSize, timeStepRate, seissolInstance.executionPlace(clusterData->size())),
      // cluster ids
      usePlasticity(usePlasticity), seissolInstance(seissolInstance), streamRuntime(4),
      globalData(&globalData), clusterData(clusterData),
      // global data
      dynRupInteriorData(dynRupInteriorData), dynRupCopyData(dynRupCopyData),
      frictionSolver(frictionSolverTemplate(ConfigVariant(Cfg()))),
      frictionSolverDevice(frictionSolverTemplateDevice(ConfigVariant(Cfg()))),
      frictionSolverCopy(frictionSolverTemplate(ConfigVariant(Cfg()))),
      frictionSolverCopyDevice(frictionSolverTemplateDevice(ConfigVariant(Cfg()))),
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

  // set timings to zero
  receiverTime = 0;

  spacetimeKernel.setGlobalData(globalData);
  timeKernel.setGlobalData(globalData);
  localKernel.setGlobalData(globalData);
  localKernel.setInitConds(&seissolInstance.getMemoryManager().getInitialConditions());
  localKernel.setGravitationalAcceleration(seissolInstance.getGravitationSetup().acceleration);
  neighborKernel.setGlobalData(globalData);
  dynamicRuptureKernel.setGlobalData(globalData);

  frictionSolver->allocateAuxiliaryMemory(globalData);
  frictionSolverDevice->allocateAuxiliaryMemory(globalData);
  frictionSolverCopy->allocateAuxiliaryMemory(globalData);
  frictionSolverCopyDevice->allocateAuxiliaryMemory(globalData);

  frictionSolver->setupLayer(*dynRupInteriorData, streamRuntime);
  frictionSolverDevice->setupLayer(*dynRupInteriorData, streamRuntime);
  frictionSolverCopy->setupLayer(*dynRupCopyData, streamRuntime);
  frictionSolverCopyDevice->setupLayer(*dynRupCopyData, streamRuntime);
  streamRuntime.wait();

  computeFlops();

  regionComputeLocalIntegration = loopStatistics->getRegion("computeLocalIntegration");
  regionComputeNeighboringIntegration = loopStatistics->getRegion("computeNeighboringIntegration");
  regionComputeDynamicRupture = loopStatistics->getRegion("computeDynamicRupture");
  regionComputePointSources = loopStatistics->getRegion("computePointSources");

  yieldCells[0] = 0;
}

template <typename Cfg>
void TimeCluster<Cfg>::setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster) {
  this->sourceCluster = std::move(sourceCluster);
}

template <typename Cfg>
void TimeCluster<Cfg>::writeReceivers() {
  SCOREP_USER_REGION("writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION)

  if (receiverCluster != nullptr) {
    receiverTime = receiverCluster->calcReceivers(
        receiverTime, ct.correctionTime, timeStepSize(), executor, streamRuntime);
  }
}

template <typename Cfg>
std::vector<NeighborCluster>* TimeCluster<Cfg>::getNeighborClusters() {
  return &neighbors;
}

template <typename Cfg>
void TimeCluster<Cfg>::computeSources() {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeSources", device::ProfilingColors::Blue);
#endif
  SCOREP_USER_REGION("computeSources", SCOREP_USER_REGION_TYPE_FUNCTION)

  // Return when point sources not initialized. This might happen if there
  // are no point sources on this rank.
  auto* pointSourceCluster = [&]() -> kernels::PointSourceCluster* {
    if (executor == Executor::Device) {
      return sourceCluster.device.get();
    } else {
      return sourceCluster.host.get();
    }
  }();

  if (pointSourceCluster != nullptr) {
    loopStatistics->begin(regionComputePointSources);
    const auto timeStepSizeLocal = timeStepSize();
    pointSourceCluster->addTimeIntegratedPointSources(
        ct.correctionTime, ct.correctionTime + timeStepSizeLocal, streamRuntime);
    loopStatistics->end(regionComputePointSources, pointSourceCluster->size(), profilingId);
  }
#ifdef ACL_DEVICE
  device.api->popLastProfilingMark();
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

  loopStatistics->begin(regionComputeDynamicRupture);

  auto* faceInformation = layerData.var<DynamicRupture::FaceInformation>();
  auto* godunovData = layerData.var<DynamicRupture::GodunovData>(Cfg());
  auto* drEnergyOutput = layerData.var<DynamicRupture::DREnergyOutputVar>(Cfg());
  auto* const* timeDerivativePlus = layerData.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  auto* const* timeDerivativeMinus = layerData.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  auto* qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(Cfg());
  auto* qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(Cfg());

  const auto timestep = timeStepSize();

  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

  const auto pointsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, timestep);
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

template <typename Cfg>
void TimeCluster<Cfg>::computeDynamicRuptureDevice(DynamicRupture::Layer& layerData) {
#ifdef ACL_DEVICE
  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeDynamicRupture);

  if (layerData.size() > 0) {
    // compute space time interpolation part

    const auto timestep = timeStepSize();

    ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, timestep);
    auto& table = layerData.getConditionalTable<inner_keys::Dr>();

    const auto [timePoints, timeWeights] =
        seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

    const auto pointsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, timestep);
    const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT(timePoints);

    streamRuntime.runGraph(
        computeGraphKey, layerData, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
          dynamicRuptureKernel.batchedSpaceTimeInterpolation(
              table, pointsCollocate.data(), streamRuntime);
        });
    device.api->popLastProfilingMark();

    auto& solver =
        &layerData == dynRupInteriorData ? frictionSolverDevice : frictionSolverCopyDevice;

    device.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    if (solver->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo<DynamicRupture::QInterpolatedPlus>(
          initializer::AllocationPlace::Host, streamRuntime.stream());
      layerData.varSynchronizeTo<DynamicRupture::QInterpolatedMinus>(
          initializer::AllocationPlace::Host, streamRuntime.stream());
      streamRuntime.wait();
      solver->evaluate(ct.correctionTime, frictionTime, timeWeights.data(), streamRuntime);
      layerData.varSynchronizeTo<DynamicRupture::FluxSolverMinus>(
          initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo<DynamicRupture::FluxSolverPlus>(
          initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo<DynamicRupture::ImposedStateMinus>(
          initializer::AllocationPlace::Device, streamRuntime.stream());
      layerData.varSynchronizeTo<DynamicRupture::ImposedStatePlus>(
          initializer::AllocationPlace::Device, streamRuntime.stream());
    } else {
      solver->evaluate(ct.correctionTime, frictionTime, timeWeights.data(), streamRuntime);
    }

    device.api->popLastProfilingMark();
  }
  loopStatistics->end(regionComputeDynamicRupture, layerData.size(), profilingId);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif
}

template <typename Cfg>
void TimeCluster<Cfg>::computeDynamicRuptureFlops(DynamicRupture::Layer& layerData,
                                                  std::uint64_t& nonZeroFlops,
                                                  std::uint64_t& hardwareFlops) {
  nonZeroFlops = 0;
  hardwareFlops = 0;

  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();

  for (std::size_t face = 0; face < layerData.size(); ++face) {
    std::uint64_t faceNonZeroFlops = 0;
    std::uint64_t faceHardwareFlops = 0;
    dynamicRuptureKernel.flopsGodunovState(
        faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::computeLocalIntegration(bool resetBuffers) {
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeLocalIntegration);

  // local integration buffer
  alignas(Alignment) real integrationBuffer[tensor::I<Cfg>::size()];

  // pointer for the call of the ADER-function
  real* bufferPointer = nullptr;

  real* const* buffers = clusterData->var<LTS::Buffers>(Cfg());
  real* const* derivatives = clusterData->var<LTS::Derivatives>(Cfg());
  CellMaterialData* materialData = clusterData->var<LTS::Material>();

  kernels::LocalTmp<Cfg> tmp(seissolInstance.getGravitationSetup().acceleration);

  const auto timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto integrationCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);

#ifdef _OPENMP
#pragma omp parallel for private(bufferPointer, integrationBuffer),                                \
    firstprivate(tmp) schedule(static)
#endif
  for (std::size_t cell = 0; cell < clusterData->size(); cell++) {
    auto data = clusterData->cellRef<Cfg>(cell);

    // We need to check, whether we can overwrite the buffer or if it is
    // needed by some other time cluster.
    // If we cannot overwrite the buffer, we compute everything in a temporary
    // local buffer and accumulate the results later in the shared buffer.
    const bool buffersProvided =
        data.template get<LTS::CellInformation>().ltsSetup.hasBuffers(); // buffers are provided
    const bool resetMyBuffers =
        buffersProvided &&
        (!data.template get<LTS::CellInformation>().ltsSetup.accumulateBuffers() ||
         resetBuffers); // they should be reset

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
    const CellBoundaryMapping<Cfg>(*boundaryMapping)[4] =
        clusterData->var<LTS::BoundaryMapping>(Cfg());
    localKernel.computeIntegral(bufferPointer,
                                data,
                                tmp,
                                &materialData[cell],
                                &boundaryMapping[cell],
                                ct.correctionTime,
                                timeStepWidth);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      auto& curFaceDisplacements = data.template get<LTS::FaceDisplacements>()[face];
      // Note: Displacement for freeSurfaceGravity is computed in Time.cpp
      if (curFaceDisplacements != nullptr &&
          data.template get<LTS::CellInformation>().faceTypes[face] !=
              FaceType::FreeSurfaceGravity) {
        kernel::addVelocity<Cfg> addVelocityKrnl;

        addVelocityKrnl.V3mTo2nFace = globalData->get<Cfg>().v3mTo2nFace;
        addVelocityKrnl.selectVelocity = init::selectVelocity<Cfg>::Values;
        addVelocityKrnl.faceDisplacement = data.template get<LTS::FaceDisplacements>()[face];
        addVelocityKrnl.I = bufferPointer;
        addVelocityKrnl.execute(face);
      }
    }

    // TODO: Integrate this step into the kernel
    // We've used a temporary buffer -> need to accumulate update in
    // shared buffer.
    if (!resetMyBuffers && buffersProvided) {
      assert(buffers[cell] != nullptr);

      for (std::size_t dof = 0; dof < tensor::I<Cfg>::size(); ++dof) {
        buffers[cell][dof] += integrationBuffer[dof];
      }
    }
  }

  loopStatistics->end(regionComputeLocalIntegration, clusterData->size(), profilingId);
}

template <typename Cfg>
void TimeCluster<Cfg>::computeLocalIntegrationDevice(bool resetBuffers) {

#ifdef ACL_DEVICE
  SCOREP_USER_REGION("computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  device.api->putProfilingMark("computeLocalIntegration", device::ProfilingColors::Yellow);

  loopStatistics->begin(regionComputeLocalIntegration);

  auto& dataTable = clusterData->getConditionalTable<inner_keys::Wp>();
  auto& materialTable = clusterData->getConditionalTable<inner_keys::Material>();
  auto& indicesTable = clusterData->getConditionalTable<inner_keys::Indices>();

  kernels::LocalTmp<Cfg> tmp(seissolInstance.getGravitationSetup().acceleration);

  const double timeStepWidth = timeStepSize();
  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
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
            dataTable, materialTable, indicesTable, timeStepWidth, streamRuntime);

        localKernel.evaluateBatchedTimeDependentBc(
            dataTable, indicesTable, *clusterData, ct.correctionTime, timeStepWidth, streamRuntime);

        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
          if (dataTable.find(key) != dataTable.end()) {
            auto& entry = dataTable[key];
            // NOTE: integrated velocities have been computed implicitly, i.e
            // it is 6th, 7the and 8th columns of integrated dofs

            kernel::gpu_addVelocity<Cfg> displacementKrnl;
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
                tensor::I<Cfg>::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime.stream());
          } else {
            device.algorithms.accumulateBatchedData(
                const_cast<const real**>(
                    (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr()),
                (entry.get(inner_keys::Wp::Id::Buffers))->getDeviceDataPtr(),
                tensor::I<Cfg>::Size,
                (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
                streamRuntime.stream());
          }
        }
      });

  loopStatistics->end(regionComputeLocalIntegration, clusterData->size(), profilingId);
  device.api->popLastProfilingMark();
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

template <typename Cfg>
void TimeCluster<Cfg>::computeNeighboringIntegration(double subTimeStart) {
  if (usePlasticity) {
    computeNeighboringIntegrationImplementation<true>(subTimeStart);
  } else {
    computeNeighboringIntegrationImplementation<false>(subTimeStart);
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::computeNeighboringIntegrationDevice(double subTimeStart) {
#ifdef ACL_DEVICE
  device.api->putProfilingMark("computeNeighboring", device::ProfilingColors::Red);
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)
  loopStatistics->begin(regionComputeNeighboringIntegration);

  const double timeStepWidth = timeStepSize();
  auto& table = clusterData->getConditionalTable<inner_keys::Wp>();

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timeStepWidth, timeStepWidth);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timeStepWidth + subTimeStart, neighborTimestep);

  seissol::kernels::TimeCommon<Cfg>::computeBatchedIntegrals(
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
        clusterData->var<LTS::Plasticity>(seissol::initializer::AllocationPlace::Device);
    auto* isAdjustableVector =
        clusterData->var<LTS::FlagScratch>(seissol::initializer::AllocationPlace::Device);
    streamRuntime.runGraph(plasticityGraphKey,
                           *clusterData,
                           [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
                             seissol::kernels::Plasticity<Cfg>::computePlasticityBatched(
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
  loopStatistics->end(regionComputeNeighboringIntegration, clusterData->size(), profilingId);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif // ACL_DEVICE
}

template <typename Cfg>
void TimeCluster<Cfg>::computeLocalIntegrationFlops() {
  auto& flopsNonZero = accFlopsNonZero[static_cast<int>(ComputePart::Local)];
  auto& flopsHardware = accFlopsHardware[static_cast<int>(ComputePart::Local)];
  flopsNonZero = 0;
  flopsHardware = 0;

  auto* cellInformation = clusterData->var<LTS::CellInformation>();
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
            GravitationalFreeSurfaceBc<Cfg>::getFlopsDisplacementFace(
                face, cellInformation[cell].faceTypes[face]);
        flopsNonZero += nonZeroFlopsDisplacement;
        flopsHardware += hardwareFlopsDisplacement;
      }
    }
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::computeNeighborIntegrationFlops() {
  auto& flopsNonZero = accFlopsNonZero[static_cast<int>(ComputePart::Neighbor)];
  auto& flopsHardware = accFlopsHardware[static_cast<int>(ComputePart::Neighbor)];
  auto& drFlopsNonZero = accFlopsNonZero[static_cast<int>(ComputePart::DRNeighbor)];
  auto& drFlopsHardware = accFlopsHardware[static_cast<int>(ComputePart::DRNeighbor)];
  flopsNonZero = 0;
  flopsHardware = 0;
  drFlopsNonZero = 0;
  drFlopsHardware = 0;

  auto* cellInformation = clusterData->var<LTS::CellInformation>();
  auto* drMapping = clusterData->var<LTS::DRMapping>(Cfg());
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

template <typename Cfg>
void TimeCluster<Cfg>::computeFlops() {
  computeLocalIntegrationFlops();
  computeNeighborIntegrationFlops();
  computeDynamicRuptureFlops(
      *dynRupInteriorData,
      accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawInterior)],
      accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawInterior)]);
  computeDynamicRuptureFlops(*dynRupCopyData,
                             accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)],
                             accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
  seissol::kernels::Plasticity<Cfg>::flopsPlasticity(
      accFlopsNonZero[static_cast<int>(ComputePart::PlasticityCheck)],
      accFlopsHardware[static_cast<int>(ComputePart::PlasticityCheck)],
      accFlopsNonZero[static_cast<int>(ComputePart::PlasticityYield)],
      accFlopsHardware[static_cast<int>(ComputePart::PlasticityYield)]);
}

template <typename Cfg>
ActResult TimeCluster<Cfg>::act() {
  actorStateStatistics->enter(state);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics->enter(state);
  return result;
}

template <typename Cfg>
void TimeCluster<Cfg>::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  if (neighborCluster.ct.maxTimeStepSize > ct.maxTimeStepSize) {
    lastSubTime = neighborCluster.ct.correctionTime;
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::handleAdvancedCorrectionTimeMessage(const NeighborCluster& /*...*/) {
  // Doesn't do anything
}

template <typename Cfg>
void TimeCluster<Cfg>::predict() {
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

  if (executor == Executor::Device) {
    computeLocalIntegrationDevice(resetBuffers);
  } else {
    computeLocalIntegration(resetBuffers);
  }

  computeSources();

  seissolInstance.flopCounter().incrementNonZeroFlopsLocal(
      accFlopsNonZero[static_cast<int>(ComputePart::Local)]);
  seissolInstance.flopCounter().incrementHardwareFlopsLocal(
      accFlopsHardware[static_cast<int>(ComputePart::Local)]);

  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                              : seissol::initializer::AllocationPlace::Device;
    clusterData->varSynchronizeTo<LTS::BuffersDerivatives>(other, streamRuntime.stream());
  }

  streamRuntime.wait();
}

template <typename Cfg>
void TimeCluster<Cfg>::handleDynamicRupture(DynamicRupture::Layer& layerData) {
  if (executor == Executor::Device) {
    computeDynamicRuptureDevice(layerData);
  } else {
    computeDynamicRupture(layerData);
  }

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                              : seissol::initializer::AllocationPlace::Device;
    layerData.varSynchronizeTo<DynamicRupture::FluxSolverMinus>(other, streamRuntime.stream());
    layerData.varSynchronizeTo<DynamicRupture::FluxSolverPlus>(other, streamRuntime.stream());
    layerData.varSynchronizeTo<DynamicRupture::ImposedStateMinus>(other, streamRuntime.stream());
    layerData.varSynchronizeTo<DynamicRupture::ImposedStatePlus>(other, streamRuntime.stream());
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::correct() {
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
  if (layerType == HaloType::Copy) {
    if (dynamicRuptureScheduler->hasDynamicRuptureFaces()) {
      handleDynamicRupture(*dynRupCopyData);
      seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(
          accFlopsNonZero[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
      seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(
          accFlopsHardware[static_cast<int>(ComputePart::DRFrictionLawCopy)]);
    }
    dynamicRuptureScheduler->setLastCorrectionStepsCopy((ct.stepsSinceStart));
  }

  if (executor == Executor::Device) {
    computeNeighboringIntegrationDevice(subTimeStart);
  } else {
    computeNeighboringIntegration(subTimeStart);
  }

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
  if (dynamicRuptureScheduler->mayComputeFaultOutput(ct.stepsSinceStart)) {
    double time = ct.correctionTime;

    // repeat the current solution for some times---to match the existing output scheme.
    // maybe replace with just writePickpointOutput(globalClusterId, time + dt, dt); some day?

    do {
      const auto oldTime = time;
      time += dynamicRuptureScheduler->getOutputTimestep();
      const auto trueTime = std::min(time, syncTime);
      const auto trueDt = trueTime - oldTime;
      faultOutputManager->writePickpointOutput(globalClusterId, trueTime, trueDt, streamRuntime);
    } while (time * (1 + 1e-8) < ct.correctionTime + ct.maxTimeStepSize);

    dynamicRuptureScheduler->setLastFaultOutput(ct.stepsSinceStart);
  }

  streamRuntime.wait();

  // TODO(Lukas) Adjust with time step rate? Relevant is maximum cluster is not on this node
  const auto nextCorrectionSteps = ct.nextCorrectionSteps();
  if (printProgress && (((nextCorrectionSteps / timeStepRate) % 100) == 0)) {
    logInfo() << "#max-updates since sync: " << nextCorrectionSteps << " @ "
              << ct.nextCorrectionTime(syncTime);
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::reset() {
  AbstractTimeCluster::reset();
  // note: redundant computation, but it needs to be done somewhere
  neighborTimestep = timeStepSize();
  for (auto& neighbor : neighbors) {
    neighborTimestep = std::max(neighbor.ct.getTimeStepSize(), neighborTimestep);
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
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

template <typename Cfg>
unsigned int TimeCluster<Cfg>::getClusterId() const {
  return clusterId;
}

template <typename Cfg>
std::size_t TimeCluster<Cfg>::layerId() const {
  return clusterData->id();
}

template <typename Cfg>
unsigned int TimeCluster<Cfg>::getGlobalClusterId() const {
  return globalClusterId;
}

template <typename Cfg>
HaloType TimeCluster<Cfg>::getLayerType() const {
  return layerType;
}

template <typename Cfg>
void TimeCluster<Cfg>::setTime(double time) {
  AbstractTimeCluster::setTime(time);
  this->receiverTime = time;
  this->lastSubTime = time;
}

template <typename Cfg>
void TimeCluster<Cfg>::finalize() {
  sourceCluster.host.reset(nullptr);
  sourceCluster.device.reset(nullptr);
  streamRuntime.dispose();

  logDebug() << "#(time steps):" << numberOfTimeSteps;
}

template <typename Cfg>
template <bool UsePlasticity>
void TimeCluster<Cfg>::computeNeighboringIntegrationImplementation(double subTimeStart) {
  const auto clusterSize = clusterData->size();
  if (clusterSize == 0) {
    return;
  }
  SCOREP_USER_REGION("computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeNeighboringIntegration);

  auto* faceNeighbors = clusterData->var<LTS::FaceNeighbors>();
  auto* drMapping = clusterData->var<LTS::DRMapping>(Cfg());
  auto* cellInformation = clusterData->var<LTS::CellInformation>();
  auto* plasticity = clusterData->var<LTS::Plasticity>(Cfg());
  auto* pstrain = clusterData->var<LTS::PStrain>(Cfg());

  // NOLINTNEXTLINE
  std::size_t numberOfTetsWithPlasticYielding = 0;

  real* timeIntegrated[4];
  real* faceNeighborsPrefetch[4];

  const auto tV = seissolInstance.getSeisSolParameters().model.tv;

  const auto timestep = timeStepSize();
  const auto oneMinusIntegratingFactor =
      seissol::kernels::Plasticity<Cfg>::computeRelaxTime(tV, timestep);

  const auto timeBasis = seissol::kernels::timeBasis<Cfg>();
  const auto timeCoeffs = timeBasis.integrate(0, timestep, timestep);
  const auto subtimeCoeffs =
      timeBasis.integrate(subTimeStart, timestep + subTimeStart, neighborTimestep);

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
               clusterData,                                                                        \
               timestep,                                                                           \
               clusterSize) reduction(+ : numberOfTetsWithPlasticYielding)
#endif
  for (std::size_t cell = 0; cell < clusterSize; cell++) {
    auto data = clusterData->cellRef<Cfg>(cell);

    real** faceNeighborsCell = reinterpret_cast<real**>(faceNeighbors[cell]);
    real** faceNeighborsCell1 = reinterpret_cast<real**>(faceNeighbors[cell + 1]);

    seissol::kernels::TimeCommon<Cfg>::computeIntegrals(
        timeKernel,
        data.template get<LTS::CellInformation>(),
        timeCoeffs.data(),
        subtimeCoeffs.data(),
        faceNeighbors[cell],
        *reinterpret_cast<real(*)[4][tensor::I<Cfg>::size()]>(
            &(globalData->get<Cfg>()
                  .integrationBufferLTS[OpenMP::threadId() * 4 *
                                        static_cast<size_t>(tensor::I<Cfg>::size())])),
        timeIntegrated);

    faceNeighborsPrefetch[0] = (cellInformation[cell].faceTypes[1] != FaceType::DynamicRupture)
                                   ? faceNeighborsCell[1]
                                   : drMapping[cell][1].godunov;
    faceNeighborsPrefetch[1] = (cellInformation[cell].faceTypes[2] != FaceType::DynamicRupture)
                                   ? faceNeighborsCell[2]
                                   : drMapping[cell][2].godunov;
    faceNeighborsPrefetch[2] = (cellInformation[cell].faceTypes[3] != FaceType::DynamicRupture)
                                   ? faceNeighborsCell[3]
                                   : drMapping[cell][3].godunov;

    // fourth face's prefetches
    if (cell + 1 < clusterSize) {
      faceNeighborsPrefetch[3] =
          (cellInformation[cell + 1].faceTypes[0] != FaceType::DynamicRupture)
              ? faceNeighborsCell1[0]
              : drMapping[cell + 1][0].godunov;
    } else {
      faceNeighborsPrefetch[3] = faceNeighborsCell[3];
    }

    neighborKernel.computeNeighborsIntegral(
        data, drMapping[cell], timeIntegrated, faceNeighborsPrefetch);

    if constexpr (UsePlasticity) {
      numberOfTetsWithPlasticYielding +=
          seissol::kernels::Plasticity<Cfg>::computePlasticity(oneMinusIntegratingFactor,
                                                               timestep,
                                                               tV,
                                                               *globalData,
                                                               &plasticity[cell],
                                                               data.template get<LTS::Dofs>(),
                                                               pstrain[cell]);
    }
#ifdef INTEGRATE_QUANTITIES
    seissolInstance.postProcessor().integrateQuantities(
        m_timeStepWidth, *clusterData, cell, dofs[cell]);
#endif // INTEGRATE_QUANTITIES
  }

  if constexpr (UsePlasticity) {
    yieldCells[0] += numberOfTetsWithPlasticYielding;
    seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
        clusterSize * accFlopsNonZero[static_cast<int>(ComputePart::PlasticityCheck)]);
    seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
        clusterSize * accFlopsHardware[static_cast<int>(ComputePart::PlasticityCheck)]);
  }

  loopStatistics->end(regionComputeNeighboringIntegration, clusterSize, profilingId);
}

template <typename Cfg>
void TimeCluster<Cfg>::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor == Executor::Host)) {
      clusterData->synchronizeTo(place, stream);
      if (layerType == HaloType::Interior) {
        dynRupInteriorData->synchronizeTo(place, stream);
      }
      if (layerType == HaloType::Copy) {
        dynRupCopyData->synchronizeTo(place, stream);
      }
    }
  }
}

template <typename Cfg>
void TimeCluster<Cfg>::finishPhase() {
  const auto cells = yieldCells[0];
  seissolInstance.flopCounter().incrementNonZeroFlopsPlasticity(
      cells * accFlopsNonZero[static_cast<int>(ComputePart::PlasticityYield)]);
  seissolInstance.flopCounter().incrementHardwareFlopsPlasticity(
      cells * accFlopsHardware[static_cast<int>(ComputePart::PlasticityYield)]);
  yieldCells[0] = 0;
}

template <typename Cfg>
std::string TimeCluster<Cfg>::description() const {
  const auto identifier = clusterData->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  return "compute-" + haloStr;
}

#define SEISSOL_CONFIGITER(cfg) template class TimeCluster<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::time_stepping

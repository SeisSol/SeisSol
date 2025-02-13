// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "DynamicRuptureCluster.h"
#include "SeisSol.h"
#include <AbstractAPI.h>
#include <Common/Executor.h>
#include <DataTypes/EncodedConstants.h>
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/Output/OutputManager.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/PreProcessorMacros.h>
#include <Initializer/Tree/Layer.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Precision.h>
#include <Monitoring/ActorStateStatistics.h>
#include <Monitoring/LoopStatistics.h>
#include <Parallel/Helper.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/ActorState.h>
#include <xdmfwriter/scorep_wrapper.h>

namespace seissol::solver::clustering::computation {

DynamicRuptureCluster::DynamicRuptureCluster(
    double maxTimeStepSize,
    long timeStepRate,
    unsigned profilingId,
    LayerType layerType,
    seissol::initializer::Layer* layer,
    seissol::initializer::DynamicRupture* descr,
    CompoundGlobalData globalData,
    seissol::dr::friction_law::FrictionSolver* frictionSolver,
    seissol::dr::friction_law::FrictionSolver* frictionSolverDevice,
    seissol::dr::output::OutputManager* faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    ActorStateStatistics* actorStateStatistics,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority)
    : FaceCluster(maxTimeStepSize,
                  timeStepRate,
                  seissolInstance.executionPlace(layer->getNumberOfCells()),
                  cpuExecutor,
                  priority),
      profilingId(profilingId), layer(layer), descr(descr), frictionSolver(frictionSolver),
      frictionSolverDevice(frictionSolverDevice), faultOutputManager(faultOutputManager),
      globalDataOnHost(globalData.onHost), globalDataOnDevice(globalData.onDevice),
      loopStatistics(loopStatistics), actorStateStatistics(actorStateStatistics),
      seissolInstance(seissolInstance), layerType(layerType) {
  dynamicRuptureKernel.setGlobalData(globalData);

  computeFlops();

  initFrictionSolverDevice();
  // TODO(David): activate here already
  // frictionSolver->copyLtsTreeToLocal(*layer, descr, 0);
  // frictionSolverDevice->copyLtsTreeToLocal(*layer, descr, 0);
  regionComputeDynamicRupture = loopStatistics->getRegion("computeDynamicRupture");
}

void DynamicRuptureCluster::handleDynamicRupture() {
#ifdef ACL_DEVICE
  if (executor == Executor::Device) {
    computeDynamicRuptureDevice();
  } else {
    computeDynamicRupture();
  }

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                              : seissol::initializer::AllocationPlace::Device;
    layer->varSynchronizeTo(descr->fluxSolverMinus, other, streamRuntime.stream());
    layer->varSynchronizeTo(descr->fluxSolverPlus, other, streamRuntime.stream());
    layer->varSynchronizeTo(descr->imposedStateMinus, other, streamRuntime.stream());
    layer->varSynchronizeTo(descr->imposedStatePlus, other, streamRuntime.stream());
  }
#else
  computeDynamicRupture();
#endif
}
void DynamicRuptureCluster::computeDynamicRupture() {
  if (layer->getNumberOfCells() == 0) {
    return;
  }

  loopStatistics->begin(regionComputeDynamicRupture);

  DRFaceInformation* faceInformation = layer->var(descr->faceInformation);
  DRGodunovData* godunovData = layer->var(descr->godunovData);
  DREnergyOutput* drEnergyOutput = layer->var(descr->drEnergyOutput);
  real** timeDerivativePlus = layer->var(descr->timeDerivativePlus);
  real** timeDerivativeMinus = layer->var(descr->timeDerivativeMinus);
  auto* qInterpolatedPlus = layer->var(descr->qInterpolatedPlus);
  auto* qInterpolatedMinus = layer->var(descr->qInterpolatedMinus);

  dynamicRuptureKernel.setTimeStepWidth(timeStepSize());
  frictionSolver->computeDeltaT(dynamicRuptureKernel.timePoints);

  const auto numberOfCells = layer->getNumberOfCells();

  streamRuntime.enqueueOmpFor(numberOfCells, [=](std::size_t face) {
    const unsigned prefetchFace = (face < numberOfCells - 1) ? face + 1 : face;
    dynamicRuptureKernel.spaceTimeInterpolation(faceInformation[face],
                                                globalDataOnHost,
                                                &godunovData[face],
                                                &drEnergyOutput[face],
                                                timeDerivativePlus[face],
                                                timeDerivativeMinus[face],
                                                qInterpolatedPlus[face],
                                                qInterpolatedMinus[face],
                                                timeDerivativePlus[prefetchFace],
                                                timeDerivativeMinus[prefetchFace]);
  });
  frictionSolver->evaluate(*layer,
                           descr,
                           ct.time[ComputeStep::Interact],
                           dynamicRuptureKernel.timeWeights,
                           streamRuntime);
  loopStatistics->end(regionComputeDynamicRupture, layer->getNumberOfCells(), profilingId);
}
void DynamicRuptureCluster::computeDynamicRuptureDevice() {
#ifdef ACL_DEVICE
  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeDynamicRupture);

  if (layer->getNumberOfCells() > 0) {
    // compute space time interpolation part

    auto& dynamicRuptureKernel = this->dynamicRuptureKernel;

    const double stepSizeWidth = timeStepSize();
    const ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, stepSizeWidth);
    auto& table = layer->getConditionalTable<inner_keys::Dr>();

    // compute constants for the kernels
    dynamicRuptureKernel.setTimeStepWidth(stepSizeWidth);
    frictionSolverDevice->computeDeltaT(dynamicRuptureKernel.timePoints);

    streamRuntime.runGraph(
        computeGraphKey, *layer, [&](seissol::parallel::runtime::StreamRuntime& streamRuntime) {
          dynamicRuptureKernel.batchedSpaceTimeInterpolation(table, streamRuntime);
        });
    device.api->popLastProfilingMark();
    if (frictionSolverDevice->allocationPlace() == initializer::AllocationPlace::Host) {
      layer->varSynchronizeTo(
          descr->qInterpolatedPlus, initializer::AllocationPlace::Host, streamRuntime.stream());
      layer->varSynchronizeTo(
          descr->qInterpolatedMinus, initializer::AllocationPlace::Host, streamRuntime.stream());
    }

    device.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);

    frictionSolverDevice->evaluate(*layer,
                                   descr,
                                   ct.time[ComputeStep::Interact],
                                   dynamicRuptureKernel.timeWeights,
                                   streamRuntime);
    device.api->popLastProfilingMark();
    if (frictionSolverDevice->allocationPlace() == initializer::AllocationPlace::Host) {
      layer->varSynchronizeTo(
          descr->fluxSolverMinus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layer->varSynchronizeTo(
          descr->fluxSolverPlus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layer->varSynchronizeTo(
          descr->imposedStateMinus, initializer::AllocationPlace::Device, streamRuntime.stream());
      layer->varSynchronizeTo(
          descr->imposedStatePlus, initializer::AllocationPlace::Device, streamRuntime.stream());
    }
  }
  loopStatistics->end(regionComputeDynamicRupture, layer->getNumberOfCells(), profilingId);
#endif
}

void DynamicRuptureCluster::computeDynamicRuptureFlops(long long& nonZeroFlops,
                                                       long long& hardwareFlops) {
  nonZeroFlops = 0;
  hardwareFlops = 0;

  DRFaceInformation* faceInformation = layer->var(descr->faceInformation);

  for (unsigned face = 0; face < layer->getNumberOfCells(); ++face) {
    long long faceNonZeroFlops = 0;
    long long faceHardwareFlops = 0;
    dynamicRuptureKernel.flopsGodunovState(
        faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

void DynamicRuptureCluster::computeFlops() {
  computeDynamicRuptureFlops(flopsNonZero, flopsHardware);
}

void DynamicRuptureCluster::runCompute(ComputeStep step) {
  if (step != ComputeStep::Interact) {
    return;
  }

  handleDynamicRupture();

  seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(flopsNonZero);
  seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(flopsHardware);

  // maybe output
  if (computePickpoint) {
    faultOutputManager->writePickpointOutput(
        ct.time[ComputeStep::Interact] + timeStepSize(), timeStepSize(), streamRuntime);
  }
}

void DynamicRuptureCluster::synchronizeTo(seissol::initializer::AllocationPlace place,
                                          void* stream) {
#ifdef ACL_DEVICE
  if ((place == initializer::AllocationPlace::Host && executor == Executor::Device) ||
      (place == initializer::AllocationPlace::Device && executor == Executor::Host)) {
    layer->synchronizeTo(place, stream);
  }
#endif
}

} // namespace seissol::solver::clustering::computation

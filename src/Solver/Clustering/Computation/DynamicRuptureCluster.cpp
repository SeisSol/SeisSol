#include "DynamicRuptureCluster.hpp"
#include "SeisSol.h"

namespace seissol::time_stepping {

DynamicRuptureCluster::DynamicRuptureCluster(
    double maxTimeStepSize,
    long timeStepRate,
    unsigned profilingId,
    seissol::initializer::Layer* layer,
    seissol::initializer::DynamicRupture* descr,
    CompoundGlobalData globalData,
    seissol::dr::friction_law::FrictionSolver* frictionSolver,
    seissol::dr::friction_law::FrictionSolver* frictionSolverDevice,
    seissol::dr::output::OutputManager* faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    ActorStateStatistics* actorStateStatistics)
    : FaceCluster(maxTimeStepSize,
                  timeStepRate,
#ifdef ACL_DEVICE
                  layer->getNumberOfCells() >= deviceHostSwitch() ? Executor::Device
                                                                  : Executor::Host
#else
                  Executor::Host
#endif
                  ),
      profilingId(profilingId), layer(layer), descr(descr), frictionSolver(frictionSolver),
      frictionSolverDevice(frictionSolverDevice), faultOutputManager(faultOutputManager),
      globalDataOnHost(globalData.onHost), globalDataOnDevice(globalData.onDevice),
      loopStatistics(loopStatistics), actorStateStatistics(actorStateStatistics),
      seissolInstance(seissolInstance) {
  dynamicRuptureKernel.setGlobalData(globalData);

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
  computeDynamicRupture(layerData);
#endif
}
void DynamicRuptureCluster::computeDynamicRupture() {
  if (layer->getNumberOfCells() == 0)
    return;

  streamRuntime.enqueueHost([=]() {
    SCOREP_USER_REGION_DEFINE(myRegionHandle)
    SCOREP_USER_REGION_BEGIN(myRegionHandle,
                             "computeDynamicRuptureSpaceTimeInterpolation",
                             SCOREP_USER_REGION_TYPE_COMMON)

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

#pragma omp parallel
    { LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation"); }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned face = 0; face < layer->getNumberOfCells(); ++face) {
      unsigned prefetchFace = (face < layer->getNumberOfCells() - 1) ? face + 1 : face;
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
    }
    SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
    {
      LIKWID_MARKER_STOP("computeDynamicRuptureSpaceTimeInterpolation");
      LIKWID_MARKER_START("computeDynamicRuptureFrictionLaw");
    }

    SCOREP_USER_REGION_BEGIN(
        myRegionHandle, "computeDynamicRuptureFrictionLaw", SCOREP_USER_REGION_TYPE_COMMON)
    frictionSolver->evaluate(*layer,
                             descr,
                             ct.time[ComputeStep::Interact],
                             dynamicRuptureKernel.timeWeights,
                             streamRuntime);
    SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
    {
      LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
    }

    loopStatistics->end(regionComputeDynamicRupture, layer->getNumberOfCells(), profilingId);
  });
}
void DynamicRuptureCluster::computeDynamicRuptureDevice() {
#ifdef ACL_DEVICE
  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics->begin(regionComputeDynamicRupture);

  if (layer->getNumberOfCells() > 0) {
    // compute space time interpolation part

    auto& dynamicRuptureKernel = this->dynamicRuptureKernel;

    const double stepSizeWidth = timeStepSize();
    ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device.api->putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    auto computeGraphKey = initializer::GraphKey(graphType, stepSizeWidth);
    auto& table = layer->getConditionalTable<inner_keys::Dr>();
    streamRuntime.enqueueHost(
        [=, &dynamicRuptureKernel]() { dynamicRuptureKernel.setTimeStepWidth(stepSizeWidth); });
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
    streamRuntime.enqueueHost([=, &dynamicRuptureKernel]() {
      frictionSolverDevice->computeDeltaT(dynamicRuptureKernel.timePoints);
    });
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
    long long faceNonZeroFlops, faceHardwareFlops;
    dynamicRuptureKernel.flopsGodunovState(
        faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

void DynamicRuptureCluster::computeFlops() {
  computeDynamicRuptureFlops(flops_nonZero, flops_hardware);
}

void DynamicRuptureCluster::runCompute(ComputeStep step) {
  if (step != ComputeStep::Interact) {
    return;
  }

  handleDynamicRupture();

  seissolInstance.flopCounter().incrementNonZeroFlopsDynamicRupture(flops_nonZero);
  seissolInstance.flopCounter().incrementHardwareFlopsDynamicRupture(flops_hardware);

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

} // namespace seissol::time_stepping

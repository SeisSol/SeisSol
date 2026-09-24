// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRuptureCluster.h"

#include "Common/Constants.h"
#include "Common/Executor.h"
#include "Common/Marker.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Numerical/Quadrature.h"
#include "SeisSol.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/FaceCluster.h"
#include "Solver/TimeStepping/StepParams.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/DeviceGraph.h"

#include <Device/AbstractAPI.h>
#endif

namespace seissol::time_stepping {

DynamicRuptureCluster::DynamicRuptureCluster(
    double maxTimeStepSize,
    long timeStepRate,
    Executor executor,
    unsigned int profilingId,
    double outputTimestep,
    CompoundGlobalData globalData,
    DynamicRupture::Layer* layerData,
    dr::friction_law::FrictionSolver* frictionSolverTemplate,
    dr::friction_law::FrictionSolver* frictionSolverTemplateDevice,
    dr::output::OutputManager* faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    ActorStateStatistics* actorStateStatistics)
    : FaceCluster(maxTimeStepSize, timeStepRate, executor), seissolInstance_(seissolInstance),
      streamRuntime_(4), globalDataOnHost_(globalData.onHost),
      globalDataOnDevice_(globalData.onDevice), layerData_(layerData),
      frictionSolver_(frictionSolverTemplate->clone()),
      frictionSolverDevice_(frictionSolverTemplateDevice->clone()),
      faultOutputManager_(faultOutputManager), outputTimestep_(outputTimestep),
      loopStatistics_(loopStatistics), actorStateStatistics_(actorStateStatistics),
      profilingId_(profilingId) {
  assert(layerData_ != nullptr);
  assert(globalDataOnHost_ != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalDataOnDevice_ != nullptr);
  }

  dynamicRuptureKernel_.setGlobalData(globalData);

  frictionSolver_->allocateAuxiliaryMemory(globalDataOnHost_);
  if constexpr (seissol::isDeviceOn()) {
    frictionSolverDevice_->allocateAuxiliaryMemory(globalDataOnDevice_);
    frictionSolverDevice_->setClock(clock_.device());
  }

  frictionSolver_->setupLayer(*layerData_, streamRuntime_);
  if constexpr (seissol::isDeviceOn()) {
    frictionSolverDevice_->setupLayer(*layerData_, streamRuntime_);
  }
  streamRuntime_.wait();

  estimate_ = computeFlops();

  regionComputeDynamicRupture_ = loopStatistics_->getRegion("computeDynamicRupture");

  perfHandle_ = seissolInstance.flopCounter().addMetric(
      getLayerType() == HaloType::Copy ? "dr-frictionlaw-copy" : "dr-frictionlaw-interior", "DR");
}

void DynamicRuptureCluster::computeDynamicRupture(const StepParams& params) {
  auto& layerData = *layerData_;
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

  const auto timestep = params.timeStepSize;

  const auto [timePoints, timeWeights] =
      seissol::quadrature::ShiftedGaussLegendre(ConvergenceOrder, 0, timestep);

  const auto pointsCollocate = seissol::kernels::timeBasis().collocate(timePoints, timestep);
  const auto frictionTime = seissol::dr::friction_law::FrictionSolver::computeDeltaT(timePoints);

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
  frictionSolver_->evaluate(params.time, frictionTime, timeWeights.data(), streamRuntime_);
  SCOREP_USER_REGION_END(myRegionHandle)
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }

  loopStatistics_->end(regionComputeDynamicRupture_, layerData.size(), profilingId_);
}

void DynamicRuptureCluster::computeDynamicRuptureDevice(
    SEISSOL_GPU_PARAM const StepParams& params) {
#ifdef ACL_DEVICE
  auto& layerData = *layerData_;

  using namespace seissol::recording;

  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  loopStatistics_->begin(regionComputeDynamicRupture_);

  if (layerData.size() > 0) {
    // compute space time interpolation part

    const auto timestep = params.timeStepSize;

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

    auto& solver = frictionSolverDevice_;

    device_.api->putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
    if (solver->allocationPlace() == initializer::AllocationPlace::Host) {
      layerData.varSynchronizeTo<DynamicRupture::QInterpolatedPlus>(
          initializer::AllocationPlace::Host, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::QInterpolatedMinus>(
          initializer::AllocationPlace::Host, streamRuntime_.stream());
      streamRuntime_.wait();
      solver->evaluate(params.time, frictionTime, timeWeights.data(), streamRuntime_);
      layerData.varSynchronizeTo<DynamicRupture::FluxSolverMinus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::FluxSolverPlus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::ImposedStateMinus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
      layerData.varSynchronizeTo<DynamicRupture::ImposedStatePlus>(
          initializer::AllocationPlace::Device, streamRuntime_.stream());
    } else {
      solver->evaluate(params.time, frictionTime, timeWeights.data(), streamRuntime_);
    }

    device_.api->popLastProfilingMark();
  }
  loopStatistics_->end(regionComputeDynamicRupture_, layerData.size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif
}

PerformanceEstimate DynamicRuptureCluster::computeFlops() {
  const auto& layerData = *layerData_;
  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();

  PerformanceEstimate estimate{};

  for (std::size_t face = 0; face < layerData.size(); ++face) {
    estimate += dynamicRuptureKernel_.metrics(faceInformation[face]);
  }

  return estimate;
}

void DynamicRuptureCluster::planPickpointOutput(const StepParams& params) {
  // repeat the current solution for some times---to match the existing output scheme.
  // maybe replace with just writePickpointOutput(layerId(), time + dt, dt); some day?
  pickpointTimes_.clear();
  double time = params.time;
  do {
    const auto oldTime = time;
    time += outputTimestep_;
    const auto trueTime = std::min(time, syncTime_);
    const auto trueDt = trueTime - oldTime;
    if (faultOutputManager_->beginPickpointStep(layerData_->id(), trueTime, trueDt)) {
      pickpointTimes_.push_back(trueTime);
    }

    // write until we've completed the current copy interval, or if we've hit a sync point
  } while (time * (1 + 1e-8) < params.time + ct_.maxTimeStepSize && time < syncTime_);
}

void DynamicRuptureCluster::writePickpointOutput(const StepParams& params) {
  const double meshDt = ct_.getTimeStepSize();
  // the friction law has just evaluated this step up to its end, and that is the state written out
  const double stateTime = params.time + params.timeStepSize;
  for (const auto time : pickpointTimes_) {
    faultOutputManager_->recordPickpointOutput(
        layerData_->id(), stateTime, time, meshDt, 0, streamRuntime_);
  }
}

StepWork DynamicRuptureCluster::prepare(ActorAction action) {
  StepWork work;
  work.hostWork = executor_ == Executor::Host || hasDifferentExecutorNeighbor() ||
                  frictionSolverDevice_->allocationPlace() == initializer::AllocationPlace::Host;
  if (action == ActorAction::Correct) {
    const auto params = stepParams();

    // without a device, the clock is up to date on the host; it has to agree with the time
    assert(isDeviceOn() || *clock_.host() == params.time);

    if (layerData_->size() > 0) {
      planPickpointOutput(params);
      work.outputs = !pickpointTimes_.empty();
      seissolInstance_.flopCounter().incrementMetric(perfHandle_, estimate_);
    }
  }
  return work;
}

void DynamicRuptureCluster::interact(const StepParams& params) {
  if (layerData_->size() == 0) {
    clock_.advance(params.timeStepSize, streamRuntime_);
    return;
  }

  if (executor_ == Executor::Device) {
    computeDynamicRuptureDevice(params);
  } else {
    computeDynamicRupture(params);
  }

  writePickpointOutput(params);

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    auto other = executor_ == Executor::Device ? seissol::initializer::AllocationPlace::Host
                                               : seissol::initializer::AllocationPlace::Device;
    layerData_->varSynchronizeTo<DynamicRupture::FluxSolverMinus>(other, streamRuntime_.stream());
    layerData_->varSynchronizeTo<DynamicRupture::FluxSolverPlus>(other, streamRuntime_.stream());
    layerData_->varSynchronizeTo<DynamicRupture::ImposedStateMinus>(other, streamRuntime_.stream());
    layerData_->varSynchronizeTo<DynamicRupture::ImposedStatePlus>(other, streamRuntime_.stream());
  }

  // the time of the cluster advances by the same step after the interaction
  clock_.advance(params.timeStepSize, streamRuntime_);

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

void* DynamicRuptureCluster::recordActionEvent() { return streamRuntime_.eventRecord(); }

void DynamicRuptureCluster::waitForEvent(SEISSOL_GPU_PARAM void* event) {
#ifdef ACL_DEVICE
  if (executor_ == Executor::Host) {
    // the host kernels run right away, so the host has to wait
    device_.api->syncEventWithHost(event);
  } else {
    streamRuntime_.eventSync(event);
  }
#endif
}

ActResult DynamicRuptureCluster::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

void DynamicRuptureCluster::finalize() {
  clock_.dispose();
  streamRuntime_.dispose();
}

void DynamicRuptureCluster::timeSet(double time) { clock_.set(time, streamRuntime_); }

void DynamicRuptureCluster::synchronizeTo(seissol::initializer::AllocationPlace place,
                                          void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor_ == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor_ == Executor::Host)) {
      layerData_->synchronizeTo(place, stream);
    }
  }
}

std::size_t DynamicRuptureCluster::layerId() const { return layerData_->id(); }

HaloType DynamicRuptureCluster::getLayerType() const { return layerData_->getIdentifier().halo; }

std::string DynamicRuptureCluster::description() const {
  const std::string haloStr = getLayerType() == HaloType::Interior ? "interior" : "copy";
  return "dr-" + haloStr;
}

} // namespace seissol::time_stepping

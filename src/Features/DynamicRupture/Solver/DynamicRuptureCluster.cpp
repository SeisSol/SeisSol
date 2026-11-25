// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRuptureCluster.h"

#include "Alignment.h"
#include "Config.h"
#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Instrumentation.h"
#include "Core/Common/Marker.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Numerical/Quadrature.h"
#include "Core/Parallel/Runtime/HostRuntime.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Parallel/Trace.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "Core/TimeStepping/Compute/FaceCluster.h"
#include "Features/DynamicRupture/Output/PickpointSampling.h"
#include "Features/DynamicRupture/Solver/DynamicRupture.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "WavePropagation/Solver/Solver.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#ifdef ACL_DEVICE
#include "Core/Parallel/DeviceGraph.h"

#include <Device/AbstractAPI.h>
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

namespace seissol::solver {

namespace {
/**
 * The time from the start of a step to a sample in it. A sample at the end of the step, within the
 * tolerance of the sampling, is taken exactly there, where the time derivatives and the friction
 * state of the step end; the sampling counts a sample within its tolerance after the end of a step
 * to that step.
 */
double timeInStep(double time, double stepStart, double stepSize) {
  if (std::abs(stepStart + stepSize - time) <=
      dr::output::PickpointSampling::Tolerance * stepSize) {
    return stepSize;
  }
  return std::clamp(time - stepStart, 0.0, stepSize);
}
} // namespace

template <typename Cfg>
DynamicRuptureCluster<Cfg>::DynamicRuptureCluster(
    double maxTimeStepSize,
    long timeStepRate,
    Executor executor,
    unsigned int profilingId,
    const std::vector<PickpointPlan>& pickpoints,
    CompoundGlobalData<Cfg> globalData,
    DynamicRupture::Layer* layerData,
    const dr::friction_law::FrictionSolverFactory& frictionSolverFactory,
    const dr::friction_law::FrictionSolverFactory& frictionSolverFactoryDevice,
    RunContext& run,
    LoopStatistics* loopStatistics,
    DeviceTimeline* deviceTimeline,
    ActorStateStatistics* actorStateStatistics)
    : FaceCluster(maxTimeStepSize,
                  timeStepRate,
                  executor,
                  layerData,
                  "dr",
                  /*faces=*/true,
                  profilingId,
                  loopStatistics,
                  deviceTimeline,
                  actorStateStatistics),
      run_(run), globalData_(globalData), layerData_(layerData),
      frictionSolver_(dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactory)),
      frictionSolverDevice_(
          dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactoryDevice)) {
  assert(layerData_ != nullptr);
  pickpoints_.reserve(pickpoints.size());
  for (const auto& plan : pickpoints) {
    pickpoints_.push_back(PickpointSlot{.output = plan.output,
                                        .sampling = plan.sampling,
                                        .times = {},
                                        .stepStartHost = 0,
                                        .stepStart = nullptr,
                                        .stepStartOnDevice = nullptr});
  }
  for (auto& slot : pickpoints_) {
    slot.stepStart = &slot.stepStartHost;
  }
  assert(globalData_.onHost != nullptr);
  if constexpr (seissol::isDeviceOn()) {
    assert(globalData_.onDevice != nullptr);
  }

  dynamicRuptureKernel_.setGlobalData(globalData);

  frictionSolver_->allocateAuxiliaryMemory(globalData_.onHost);
  if constexpr (seissol::isDeviceOn()) {
    frictionSolverDevice_->allocateAuxiliaryMemory(globalData_.onDevice);
    frictionSolverDevice_->setClock(clock_.device());
  }

  frictionSolver_->setupLayer(*layerData_, streamRuntime_);
  if constexpr (seissol::isDeviceOn()) {
    frictionSolverDevice_->setupLayer(*layerData_, streamRuntime_);
  }
  streamRuntime_.wait();

  estimate_ = computeFlops();

  regionComputeDynamicRupture_ = loopStatistics_->getRegion("computeDynamicRupture");
  traceDynamicRupture_ = monitoring::Trace::name("computeDynamicRupture");

  perfHandle_ = run_.flopCounter().addMetric(
      getLayerType() == HaloType::Copy ? "dr-frictionlaw-copy" : "dr-frictionlaw-interior", "DR");
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::computeDynamicRupture(const StepParams& params) {
  auto& layerData = *layerData_;
  if (layerData.size() == 0) {
    return;
  }
  SCOREP_USER_REGION_DEFINE(myRegionHandle)
  SCOREP_USER_REGION_BEGIN(
      myRegionHandle, "computeDynamicRuptureSpaceTimeInterpolation", SCOREP_USER_REGION_TYPE_COMMON)

  const auto loopBegin = LoopStatistics::begin();

  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();
  const DRGodunovData<Cfg>* godunovData = layerData.var<DynamicRupture::GodunovData>(Cfg());
  real* const* timeDerivativePlus = layerData.var<DynamicRupture::TimeDerivativePlus>(Cfg());
  real* const* timeDerivativeMinus = layerData.var<DynamicRupture::TimeDerivativeMinus>(Cfg());
  auto* qInterpolatedPlus = layerData.var<DynamicRupture::QInterpolatedPlus>(Cfg());
  auto* qInterpolatedMinus = layerData.var<DynamicRupture::QInterpolatedMinus>(Cfg());

  const auto timestep = params.timeStepSize;

  const auto [timePoints, timeWeights] =
      seissol::quadrature::shiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

  const auto pointsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, timestep);
  const auto frictionTime =
      seissol::dr::friction_law::FrictionSolver::computeDeltaT<Cfg>(timePoints);

  // the markers of LIKWID count per thread; without LIKWID, the regions would cost a fork and a
  // join each, for nothing
#ifdef LIKWID_PERFMON
#pragma omp parallel
  {
    LIKWID_MARKER_START("computeDynamicRuptureSpaceTimeInterpolation");
  }
#endif

  const parallel::runtime::HostRuntime::LoopName loopName(traceDynamicRupture_);
  parallel::runtime::HostRuntime::forRange(
      layerData.size(), [&](std::size_t begin, std::size_t end) {
        for (std::size_t face = begin; face < end; ++face) {
          // the next face of the range comes next on this thread
          const std::size_t prefetchFace = (face + 1 < end) ? face + 1 : face;
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
      });
  SCOREP_USER_REGION_END(myRegionHandle)
#ifdef LIKWID_PERFMON
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("computeDynamicRuptureSpaceTimeInterpolation");
    LIKWID_MARKER_START("computeDynamicRuptureFrictionLaw");
  }
#endif

  SCOREP_USER_REGION_BEGIN(
      myRegionHandle, "computeDynamicRuptureFrictionLaw", SCOREP_USER_REGION_TYPE_COMMON)
  frictionSolver_->evaluate(params.time, frictionTime, timeWeights.data(), streamRuntime_);
  SCOREP_USER_REGION_END(myRegionHandle)
#ifdef LIKWID_PERFMON
#pragma omp parallel
  {
    LIKWID_MARKER_STOP("computeDynamicRuptureFrictionLaw");
  }
#endif

  loopStatistics_->end(regionComputeDynamicRupture_, layerData.size(), profilingId_, loopBegin);
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::computeDynamicRuptureDevice(
    SEISSOL_GPU_PARAM const StepParams& params) {
#ifdef ACL_DEVICE
  auto& layerData = *layerData_;

  using namespace seissol::recording;

  SCOREP_USER_REGION("computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION)

  const auto loop = deviceLoops_.begin(streamRuntime_.stream());

  if (layerData.size() > 0) {
    // compute space time interpolation part

    const auto timestep = params.timeStepSize;

    const ComputeGraphType graphType = ComputeGraphType::DynamicRuptureInterface;
    device_.api().putProfilingMark("computeDrInterfaces", device::ProfilingColors::Cyan);
    const auto computeGraphKey = initializer::GraphKey(graphType, timestep);
    auto& table = layerData.getConditionalTable();

    const auto [timePoints, timeWeights] =
        seissol::quadrature::shiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

    const auto pointsCollocate = seissol::kernels::timeBasis<Cfg>().collocate(timePoints, timestep);
    const auto frictionTime =
        seissol::dr::friction_law::FrictionSolver::computeDeltaT<Cfg>(timePoints);

    auto& solver = frictionSolverDevice_;
    // a friction law on the device takes the time of its step from the clock of the cluster, and
    // otherwise only what the step size determines; so it goes into the graph of the interpolation
    const bool frictionOnDevice = solver->allocationPlace() != initializer::AllocationPlace::Host;
    // (a lambda cannot capture the structured binding before C++20)
    const double* frictionWeights = timeWeights.data();

    streamRuntime_.runGraph(
        computeGraphKey,
        layerData,
        [&](seissol::parallel::runtime::StreamRuntime& /*streamRuntime*/) {
          dynamicRuptureKernel_.batchedSpaceTimeInterpolation(
              table, pointsCollocate.data(), streamRuntime_);
          if (frictionOnDevice) {
            solver->evaluate(params.time, frictionTime, frictionWeights, streamRuntime_);
          }
        },
        isRecurringTimestep(timestep));
    device_.api().popLastProfilingMark();

    if (!frictionOnDevice) {
      device_.api().putProfilingMark("evaluateFriction", device::ProfilingColors::Lime);
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
      device_.api().popLastProfilingMark();
    }
  }
  deviceLoops_.end(
      loop, streamRuntime_.stream(), regionComputeDynamicRupture_, layerData.size(), profilingId_);
#else
  logError() << "The GPU kernels are disabled in this version of SeisSol.";
#endif
}

template <typename Cfg>
PerformanceEstimate DynamicRuptureCluster<Cfg>::computeFlops() {
  const auto& layerData = *layerData_;
  const DRFaceInformation* faceInformation = layerData.var<DynamicRupture::FaceInformation>();

  PerformanceEstimate estimate{};

  for (std::size_t face = 0; face < layerData.size(); ++face) {
    estimate += dynamicRuptureKernel_.metrics(faceInformation[face]);
  }

  return estimate;
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::planPickpointOutput(const StepParams& params) {
  for (auto& slot : pickpoints_) {
    slot.times.clear();
    if (slot.output->hasReceivers(layerData_->id())) {
      slot.times = slot.sampling.times(params.time, params.timeStepSize);
    }
  }
}

template <typename Cfg>
std::vector<std::vector<double>> DynamicRuptureCluster<Cfg>::plannedPickpointTimes() const {
  std::vector<std::vector<double>> times;
  times.reserve(pickpoints_.size());
  for (const auto& slot : pickpoints_) {
    times.push_back(slot.times);
  }
  return times;
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::setRunTimeOutputs(bool runTimeOutputs) {
  ComputeCluster::setRunTimeOutputs(runTimeOutputs);
#ifdef ACL_DEVICE
  for (auto& slot : pickpoints_) {
    if (runTimeOutputs_ && slot.stepStart == &slot.stepStartHost) {
      slot.stepStart = static_cast<double*>(device_.api().allocPinnedMem(sizeof(double)));
      slot.stepStartOnDevice = static_cast<double*>(device_.api().devicePointer(slot.stepStart));
    }
  }
#endif
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::recordPickpointsNow(std::size_t slot,
                                                     double stepTime,
                                                     bool lastStep) {
  // the size of the step that starts at stepTime, as ClusterTimes::timeStepSize() decides it
  const double stepSize =
      lastStep ? syncTime_ - stepTime : std::min(syncTime_ - stepTime, ct_.maxTimeStepSize);
  const auto layerId = layerData_->id();
  const auto& pickpoints = pickpoints_[slot];
  for (const auto time : pickpoints.sampling.times(stepTime, stepSize)) {
    pickpoints.output->evaluate(layerId, time, stepSize, timeInStep(time, stepTime, stepSize));
  }
  pickpoints.output->finishStepNow(layerId);
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::writePickpointOutput(const StepParams& params,
                                                      const std::vector<std::vector<double>>& times,
                                                      bool lastStep) {
  const auto layerId = layerData_->id();
  for (std::size_t slot = 0; slot < pickpoints_.size(); ++slot) {
    const auto& pickpoints = pickpoints_[slot];
    auto* output = pickpoints.output;
    if (!output->hasReceivers(layerId)) {
      continue;
    }

    if (runTimeOutputs_) {
      if (output->postAtRunTime(
              layerId, clock_.device(), streamRuntime_, [this, slot, lastStep](double stepStart) {
                recordPickpointsNow(slot, stepStart, lastStep);
              })) {
        continue;
      }
      auto& callRuntime = output->gatherData(layerId,
                                             clock_.device(),
                                             pickpoints.stepStart,
                                             pickpoints.stepStartOnDevice,
                                             params.time,
                                             streamRuntime_);
      callRuntime.enqueueHost([this, slot, lastStep] {
        recordPickpointsNow(slot, *pickpoints_[slot].stepStart, lastStep);
      });
      continue;
    }

    for (const auto time : times[slot]) {
      output->record(layerId,
                     time,
                     params.timeStepSize,
                     timeInStep(time, params.time, params.timeStepSize),
                     streamRuntime_);
    }
    output->finishStep(layerId, !times[slot].empty(), streamRuntime_);
  }
}

template <typename Cfg>
std::vector<parallel::runtime::StreamRuntime*> DynamicRuptureCluster<Cfg>::runtimes() {
  auto runtimes = ComputeCluster::runtimes();
  for (const auto& slot : pickpoints_) {
    if (slot.output->hasReceivers(layerData_->id())) {
      auto* runtime = slot.output->runtime(layerData_->id());
      if (runtime != nullptr) {
        runtimes.push_back(runtime);
      }
    }
  }
  return runtimes;
}

template <typename Cfg>
bool DynamicRuptureCluster<Cfg>::hostWork() const {
  return ComputeCluster::hostWork() ||
         frictionSolverDevice_->allocationPlace() == initializer::AllocationPlace::Host;
}

template <typename Cfg>
bool DynamicRuptureCluster<Cfg>::outputsAhead(long steps) const {
  if (layerData_->size() == 0 || runTimeOutputs_) {
    return false;
  }
  // the same sample times as the host parts of the next steps plan
  for (const auto& slot : pickpoints_) {
    if (!slot.output->hasReceivers(layerData_->id())) {
      continue;
    }
    double stepTime = ct_.correctionTime;
    for (long step = 0; step < steps; ++step) {
      const double stepSize =
          ct_.timeStepSize(stepTime, ct_.stepsSinceLastSync + step * ct_.timeStepRate, syncTime_);
      if (slot.sampling.any(stepTime, stepSize)) {
        return true;
      }
      stepTime += stepSize;
    }
  }
  return false;
}

template <typename Cfg>
StepWork DynamicRuptureCluster<Cfg>::prepare(ActorAction action) {
  StepWork work;
  work.hostWork = hostWork();
  if (action == ActorAction::Correct) {
    const auto params = stepParams();

    // without a device, the clock is up to date on the host; it has to agree with the time
    assert(isDeviceOn() || *clock_.host() == params.time);

    if (layerData_->size() > 0) {
      if (!runTimeOutputs_) {
        planPickpointOutput(params);
        work.outputs =
            std::ranges::any_of(pickpoints_, [](const auto& slot) { return !slot.times.empty(); });
      }
      run_.flopCounter().incrementMetric(perfHandle_, estimate_);
    }
  }
  return work;
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::interact(const StepParams& params) {
  deviceLoops_.collect();
  if (layerData_->size() == 0) {
    clock_.advance(params.timeStepSize, streamRuntime_);
    return;
  }

  // a recording that gets replayed ends before the synchronization point, so the last step
  // before it is never replayed
  const bool lastStep = ct_.isLastStep(ct_.stepsSinceLastSync);

  if (executor_ == Executor::Device) {
    computeDynamicRuptureDevice(params);
    writePickpointOutput(params, plannedPickpointTimes(), lastStep);
  } else {
    // the faces read what the adjacent cells share, and change what these read for their
    // correction; the work takes the sample times of its step along, as the decisions about the
    // next steps may come before it runs
    hostRuntime_.action({/*readsNeighbors=*/true, /*writesShared=*/true},
                        [this, params, times = plannedPickpointTimes(), lastStep] {
                          computeDynamicRupture(params);
                          writePickpointOutput(params, times, lastStep);
                        });
  }

  // TODO(David): restrict to copy/interior of same cluster type
  if (hasDifferentExecutorNeighbor()) {
    const auto other = executor_ == Executor::Device
                           ? seissol::initializer::AllocationPlace::Host
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

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::finalize() {
#ifdef ACL_DEVICE
  for (auto& slot : pickpoints_) {
    if (slot.stepStart != &slot.stepStartHost) {
      device_.api().freePinnedMem(slot.stepStart);
      slot.stepStart = &slot.stepStartHost;
      slot.stepStartOnDevice = nullptr;
    }
  }
#endif
  ComputeCluster::finalize();
}

#define SEISSOL_INSTANTIATE(Cfg) template class DynamicRuptureCluster<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::solver

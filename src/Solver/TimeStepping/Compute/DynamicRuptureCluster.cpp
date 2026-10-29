// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRuptureCluster.h"

#include "Common/ConfigRegistry.h"
#include "Common/Executor.h"
#include "Common/Marker.h"
#include "Config.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "DynamicRupture/Output/PickpointSampling.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/DynamicRupture.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/DeviceTimeline.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Monitoring/Trace.h"
#include "Numerical/Quadrature.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "SeisSol.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Actor/StepParams.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Initializer/DeviceGraph.h"

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
    const dr::output::PickpointSampling& pickpointSampling,
    CompoundGlobalData<Cfg> globalData,
    DynamicRupture::Layer* layerData,
    const dr::friction_law::FrictionSolverFactory& frictionSolverFactory,
    const dr::friction_law::FrictionSolverFactory& frictionSolverFactoryDevice,
    dr::output::OutputManager* faultOutputManager,
    seissol::SeisSol& seissolInstance,
    LoopStatistics* loopStatistics,
    DeviceTimeline* deviceTimeline,
    ActorStateStatistics* actorStateStatistics)
    : DynamicRuptureClusterInterface(maxTimeStepSize, timeStepRate, executor),
      seissolInstance_(seissolInstance),
      streamRuntime_(
          4,
          parallel::runtime::ownQueuePriorityOf(layerData->getIdentifier(), /*faces=*/true),
          parallel::runtime::ownQueueOf(layerData->getIdentifier(), /*faces=*/true)),
      globalData_(globalData), layerData_(layerData),
      frictionSolver_(dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactory)),
      frictionSolverDevice_(dr::friction_law::makeFrictionSolver<Cfg>(frictionSolverFactoryDevice)),
      faultOutputManager_(faultOutputManager), pickpointSampling_(pickpointSampling),
      loopStatistics_(loopStatistics),
      deviceLoops_(loopStatistics, deviceTimeline, DeviceTimeline::Kind::Loops),
      actorStateStatistics_(actorStateStatistics), profilingId_(profilingId) {
  assert(layerData_ != nullptr);
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
  // (named once the trace records, after the cluster is complete)
  deviceLoops_.setTrace([this] { return description(); });

  perfHandle_ = seissolInstance.flopCounter().addMetric(
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
      seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

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
    auto computeGraphKey = initializer::GraphKey(graphType, timestep);
    auto& table = layerData.getConditionalTable<inner_keys::Dr>();

    const auto [timePoints, timeWeights] =
        seissol::quadrature::ShiftedGaussLegendre(Cfg::ConvergenceOrder, 0, timestep);

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
  pickpointTimes_.clear();
  if (faultOutputManager_->hasPickpoints(layerData_->id())) {
    pickpointTimes_ = pickpointSampling_.times(params.time, params.timeStepSize);
  }
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::setRunTimeOutputs(bool runTimeOutputs) {
  runTimeOutputs_ = runTimeOutputs;
#ifdef ACL_DEVICE
  if (runTimeOutputs_ && pickpointStepStart_ == &pickpointStepStartHost_) {
    pickpointStepStart_ = static_cast<double*>(device_.api().allocPinnedMem(sizeof(double)));
    pickpointStepStartOnDevice_ =
        static_cast<double*>(device_.api().devicePointer(pickpointStepStart_));
  }
#endif
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::recordPickpointsNow(double stepTime, bool lastStep) {
  // the size of the step that starts at stepTime, as ClusterTimes::timeStepSize() decides it
  const double stepSize =
      lastStep ? syncTime_ - stepTime : std::min(syncTime_ - stepTime, ct_.maxTimeStepSize);
  const auto layerId = layerData_->id();
  for (const auto time : pickpointSampling_.times(stepTime, stepSize)) {
    faultOutputManager_->evaluatePickpointOutput(
        layerId, time, stepSize, timeInStep(time, stepTime, stepSize));
  }
  faultOutputManager_->finishPickpointStepNow(layerId);
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::writePickpointOutput(const StepParams& params,
                                                      const std::vector<double>& times,
                                                      bool lastStep) {
  if (runTimeOutputs_) {
    const auto layerId = layerData_->id();
    if (faultOutputManager_->hasPickpoints(layerId)) {
      if (faultOutputManager_->postPickpointsAtRunTime(
              layerId, clock_.device(), streamRuntime_, [this, lastStep](double stepStart) {
                recordPickpointsNow(stepStart, lastStep);
              })) {
        return;
      }
      auto& callRuntime = faultOutputManager_->gatherPickpointData(layerId,
                                                                   clock_.device(),
                                                                   pickpointStepStart_,
                                                                   pickpointStepStartOnDevice_,
                                                                   params.time,
                                                                   streamRuntime_);
      callRuntime.enqueueHost(
          [this, lastStep]() { recordPickpointsNow(*pickpointStepStart_, lastStep); });
    }
    return;
  }

  if (!faultOutputManager_->hasPickpoints(layerData_->id())) {
    return;
  }
  for (const auto time : times) {
    faultOutputManager_->recordPickpointOutput(layerData_->id(),
                                               time,
                                               params.timeStepSize,
                                               timeInStep(time, params.time, params.timeStepSize),
                                               streamRuntime_);
  }
  faultOutputManager_->finishPickpointStep(layerData_->id(), !times.empty(), streamRuntime_);
}

template <typename Cfg>
std::vector<parallel::runtime::StreamRuntime*> DynamicRuptureCluster<Cfg>::runtimes() {
  std::vector<parallel::runtime::StreamRuntime*> runtimes{&streamRuntime_};
  if (faultOutputManager_ != nullptr && faultOutputManager_->hasPickpoints(layerData_->id())) {
    auto* runtime = faultOutputManager_->pickpointRuntime(layerData_->id());
    if (runtime != nullptr) {
      runtimes.push_back(runtime);
    }
  }
  return runtimes;
}

template <typename Cfg>
bool DynamicRuptureCluster<Cfg>::hostWork() const {
  return executor_ == Executor::Host || hasDifferentExecutorNeighbor() ||
         frictionSolverDevice_->allocationPlace() == initializer::AllocationPlace::Host;
}

template <typename Cfg>
bool DynamicRuptureCluster<Cfg>::outputsAhead(long steps) const {
  if (layerData_->size() == 0 || runTimeOutputs_ ||
      !faultOutputManager_->hasPickpoints(layerData_->id())) {
    return false;
  }
  // the same sample times as the host parts of the next steps plan
  double stepTime = ct_.correctionTime;
  for (long step = 0; step < steps; ++step) {
    const double stepSize =
        ct_.timeStepSize(stepTime, ct_.stepsSinceLastSync + step * ct_.timeStepRate, syncTime_);
    if (pickpointSampling_.any(stepTime, stepSize)) {
      return true;
    }
    stepTime += stepSize;
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
        work.outputs = !pickpointTimes_.empty();
      }
      seissolInstance_.flopCounter().incrementMetric(perfHandle_, estimate_);
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
    writePickpointOutput(params, pickpointTimes_, lastStep);
  } else {
    // the faces read what the adjacent cells share, and change what these read for their
    // correction; the work takes the sample times of its step along, as the decisions about the
    // next steps may come before it runs
    hostRuntime_.action({/*readsNeighbors=*/true, /*writesShared=*/true},
                        [this, params, times = pickpointTimes_, lastStep]() {
                          computeDynamicRupture(params);
                          writePickpointOutput(params, times, lastStep);
                        });
  }

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

template <typename Cfg>
ActorEvent DynamicRuptureCluster<Cfg>::recordActionEvent() {
  return streamRuntime_.eventRecord();
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::waitForEvent(SEISSOL_GPU_PARAM const ActorEvent& event) {
#ifdef ACL_DEVICE
  if (executor_ == Executor::Host) {
    // the host kernels run right away, so the host has to wait (they never get recorded)
    device_.api().syncEventWithHost(event.get());
  } else {
    streamRuntime_.eventSync(event);
  }
#endif
}

template <typename Cfg>
ActResult DynamicRuptureCluster<Cfg>::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::finishPhase() {
  // the device is past all loops of the epoch
  deviceLoops_.collectAll();
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::finalize() {
  deviceLoops_.release();
#ifdef ACL_DEVICE
  if (pickpointStepStart_ != &pickpointStepStartHost_) {
    device_.api().freePinnedMem(pickpointStepStart_);
    pickpointStepStart_ = &pickpointStepStartHost_;
    pickpointStepStartOnDevice_ = nullptr;
  }
#endif
  AbstractTimeCluster::finalize();
  clock_.dispose();
  streamRuntime_.dispose();
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::timeSet(double time) {
  clock_.set(time, streamRuntime_);
}

template <typename Cfg>
void DynamicRuptureCluster<Cfg>::synchronizeTo(seissol::initializer::AllocationPlace place,
                                               void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor_ == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor_ == Executor::Host)) {
      layerData_->synchronizeTo(place, stream);
    }
  }
}

template <typename Cfg>
std::size_t DynamicRuptureCluster<Cfg>::layerId() const {
  return layerData_->id();
}

template <typename Cfg>
HaloType DynamicRuptureCluster<Cfg>::getLayerType() const {
  return layerData_->getIdentifier().halo;
}

template <typename Cfg>
std::string DynamicRuptureCluster<Cfg>::description() const {
  const auto identifier = layerData_->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  // the time cluster, and the configuration where there can be several
  return "dr-" + haloStr + "-" + std::to_string(identifier.lts) +
         (builtConfigCount() > 1 ? "-config" + std::to_string(identifier.config) : std::string());
}

#define SEISSOL_INSTANTIATE(Cfg) template class DynamicRuptureCluster<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::solver

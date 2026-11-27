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
#include "Core/Common/Executor.h"
#include "Core/Common/Instrumentation.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ClusterStage.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <utility>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/AbstractAPI.h>
#endif

namespace seissol::solver {

CellCluster::CellCluster(unsigned int clusterId,
                         unsigned int globalClusterId,
                         unsigned int profilingId,
                         double maxTimeStepSize,
                         long timeStepRate,
                         Executor executor,
                         bool printProgress,
                         LTS::Layer* clusterData,
                         LoopStatistics* loopStatistics,
                         DeviceTimeline* deviceTimeline,
                         ActorStateStatistics* actorStateStatistics)
    : ComputeCluster(maxTimeStepSize,
                     timeStepRate,
                     executor,
                     clusterData,
                     "compute",
                     /*faces=*/false,
                     profilingId,
                     loopStatistics,
                     deviceTimeline,
                     actorStateStatistics),
      clusterData_(clusterData), printProgress_(printProgress), clusterId_(clusterId),
      globalClusterId_(globalClusterId) {
  assert(clusterData_ != nullptr);
}

void CellCluster::addStage(std::unique_ptr<ClusterStage> stage) {
  stage->setClock(clock_.device());
  stageRegions_.push_back(loopStatistics_->getRegion(stage->name()));
  stages_.push_back(std::move(stage));
}

void CellCluster::runStages(StagePhase phase, const StepParams& params) {
  SCOREP_USER_REGION("runStages", SCOREP_USER_REGION_TYPE_FUNCTION)
  for (std::size_t i = 0; i < stages_.size(); ++i) {
    auto& stage = *stages_[i];
    if (stage.phase() != phase || !stage.hasWork(executor_)) {
      continue;
    }
#ifdef ACL_DEVICE
    device_.api().putProfilingMark(stage.name(), device::ProfilingColors::Blue);
#endif
    if (executor_ == Executor::Device) {
      const auto loop = deviceLoops_.begin(streamRuntime_.stream());
      stage.run(params, executor_, streamRuntime_);
      deviceLoops_.end(
          loop, streamRuntime_.stream(), stageRegions_[i], stage.size(executor_), profilingId_);
    } else {
      const auto loopBegin = LoopStatistics::begin();
      stage.run(params, executor_, streamRuntime_);
      loopStatistics_->end(stageRegions_[i], stage.size(executor_), profilingId_, loopBegin);
    }
#ifdef ACL_DEVICE
    device_.api().popLastProfilingMark();
#endif
  }
}

StepWork CellCluster::prepare(ActorAction action) {
  const auto params = stepParams();
  StepWork work;
  work.hostWork = hostWork();

  if (action == ActorAction::Predict) {
    // without a device, the clock is up to date on the host; it has to agree with the time
    assert(isDeviceOn() || *clock_.host() == params.time);

    if (clusterData_->size() > 0) {
      work.outputs = planPrediction(params);
    }
  }

  if (action == ActorAction::Correct) {
    planCorrection(params);

    if (printProgress_) {
      const auto nextCorrectionSteps = ct_.nextCorrectionSteps();
      if (((nextCorrectionSteps / timeStepRate_) % 100) == 0) {
        logInfo() << "Max cluster / LTS cycle updates since sync: " << nextCorrectionSteps
                  << " at time " << ct_.nextCorrectionTime(syncTime_);
      }
    }
  }
  return work;
}

void CellCluster::predict() {
  assert(state_ == ActorState::Corrected);
  const auto params = stepParams();
  deviceLoops_.collect();

  if (clusterData_->size() == 0) {
    return;
  }

  if (executor_ == Executor::Device) {
    devicePrediction(params);
    runStages(StagePhase::AfterPrediction, params);
  } else {
    // the prediction changes what the neighbors read; the work takes what was decided about its
    // step along, as the decisions about the next steps may come before it runs
    hostRuntime_.action({.readsNeighbors = false, .writesShared = true},
                        [this, params, prediction = hostPrediction(params)] {
                          prediction();
                          runStages(StagePhase::AfterPrediction, params);
                        });
  }

  if (hasDifferentExecutorNeighbor()) {
    const auto other = executor_ == Executor::Device
                           ? seissol::initializer::AllocationPlace::Host
                           : seissol::initializer::AllocationPlace::Device;
    synchronizeShared(other);
  }

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

void CellCluster::correct() {
  assert(state_ == ActorState::Predicted);
  const auto params = stepParams();
  deviceLoops_.collect();

  if (executor_ == Executor::Device) {
    deviceCorrection(params);
  } else {
    // the correction reads what the neighbors share, and changes the DOFs only
    hostRuntime_.action({.readsNeighbors = true, .writesShared = false},
                        [this, params] { hostCorrection(params); });
  }

  // the time of the cluster advances by the same step after the correction
  clock_.advance(params.timeStepSize, streamRuntime_);

  if (!concurrent()) {
    streamRuntime_.wait();
  }
}

void CellCluster::predictOnce(const StepParams& params) {
  // (the device timer would keep the events of the loops otherwise)
  deviceLoops_.collect();
  if (executor_ == Executor::Device) {
    devicePrediction(params);
  } else {
    hostPrediction(params)();
  }
}

void CellCluster::correctOnce(const StepParams& params) {
  deviceLoops_.collect();
  if (executor_ == Executor::Device) {
    deviceCorrection(params);
  } else {
    hostCorrection(params);
  }
}

void CellCluster::finalize() {
  stages_.clear();
  ComputeCluster::finalize();

  logDebug() << "#(time steps):" << numberOfTimeSteps_;
}

} // namespace seissol::solver

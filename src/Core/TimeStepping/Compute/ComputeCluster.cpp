// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ComputeCluster.h"

#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Marker.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/DeviceTimer.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace seissol::solver {

ComputeCluster::ComputeCluster(double maxTimeStepSize,
                               long timeStepRate,
                               Executor executor,
                               initializer::Layer* layer,
                               std::string kind,
                               bool faces,
                               unsigned int profilingId,
                               LoopStatistics* loopStatistics,
                               DeviceTimeline* deviceTimeline,
                               ActorStateStatistics* actorStateStatistics)
    : AbstractTimeCluster(maxTimeStepSize, timeStepRate, executor), layer_(layer),
      streamRuntime_(4,
                     parallel::runtime::hsaQueuePriorityOf(layer->getIdentifier(), faces),
                     parallel::runtime::hsaQueueOf(layer->getIdentifier(), faces)),
      loopStatistics_(loopStatistics),
      deviceLoops_(loopStatistics, deviceTimeline, DeviceTimeline::Kind::Loops),
      actorStateStatistics_(actorStateStatistics), profilingId_(profilingId),
      kind_(std::move(kind)) {
  // (named once the trace records, after the cluster is complete)
  deviceLoops_.setTrace([this] { return description(); });
}

ActResult ComputeCluster::act() {
  actorStateStatistics_->enter(state_);
  const auto result = AbstractTimeCluster::act();
  actorStateStatistics_->enter(state_);
  return result;
}

void ComputeCluster::finalize() {
  deviceLoops_.release();
  AbstractTimeCluster::finalize();
  clock_.dispose();
  streamRuntime_.dispose();
}

void ComputeCluster::finishPhase() {
  // the device is past all loops of the epoch
  deviceLoops_.collectAll();
}

void ComputeCluster::synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) {
  if constexpr (isDeviceOn()) {
    if ((place == initializer::AllocationPlace::Host && executor_ == Executor::Device) ||
        (place == initializer::AllocationPlace::Device && executor_ == Executor::Host)) {
      layer_->synchronizeTo(place, stream);
    }
  }
}

void ComputeCluster::setRunTimeOutputs(bool runTimeOutputs) { runTimeOutputs_ = runTimeOutputs; }

const double* ComputeCluster::deviceTime() const { return isDeviceOn() ? clock_.host() : nullptr; }

bool ComputeCluster::hostWork() const {
  return executor_ == Executor::Host || hasDifferentExecutorNeighbor();
}

std::vector<parallel::runtime::StreamRuntime*> ComputeCluster::runtimes() {
  return {&streamRuntime_};
}

std::string ComputeCluster::description() const {
  const auto identifier = layer_->getIdentifier();
  const std::string haloStr = identifier.halo == HaloType::Interior ? "interior" : "copy";
  // the time cluster, and the configuration where there can be several
  return kind_ + "-" + haloStr + "-" + std::to_string(identifier.lts) +
         (builtConfigCount() > 1 ? "-config" + std::to_string(identifier.config) : std::string());
}

std::size_t ComputeCluster::layerId() const { return layer_->id(); }

HaloType ComputeCluster::getLayerType() const { return layer_->getIdentifier().halo; }

ActorEvent ComputeCluster::recordActionEvent() { return streamRuntime_.eventRecord(); }

void ComputeCluster::waitForEvent(SEISSOL_GPU_PARAM const ActorEvent& event) {
#ifdef ACL_DEVICE
  if (executor_ == Executor::Host) {
    // the host kernels run right away, so the host has to wait (they never get recorded)
    device_.api().syncEventWithHost(event.get());
  } else {
    streamRuntime_.eventSync(event);
  }
#endif
}

void ComputeCluster::timeSet(double time) { clock_.set(time, streamRuntime_); }

} // namespace seissol::solver

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_COMPUTECLUSTER_H_
#define SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_COMPUTECLUSTER_H_

#include "Core/Common/BasicTypedefs.h"
#include "Core/Common/Executor.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Monitoring/ActorStateStatistics.h"
#include "Core/Monitoring/DeviceTimer.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Compute/ClusterClock.h"

#include <cstddef>
#include <string>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

/**
 * The frame of a cluster that computes on a layer of a storage: its stream and clock, the events of
 * its actions, the statistics of its loops and of its actor states, and its name. The cell and the
 * face clusters derive from it, and do the work of the steps.
 */
class ComputeCluster : public AbstractTimeCluster {
  public:
  ~ComputeCluster() override = default;

  ComputeCluster(const ComputeCluster&) = delete;
  ComputeCluster(ComputeCluster&&) = delete;
  auto operator=(const ComputeCluster&) = delete;
  auto operator=(ComputeCluster&&) = delete;

  ActResult act() override;
  void finalize() override;
  void finishPhase() override;
  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;
  void setRunTimeOutputs(bool runTimeOutputs) override;

  [[nodiscard]] const double* deviceTime() const override;
  [[nodiscard]] bool hostWork() const override;
  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() override;
  [[nodiscard]] std::string description() const override;

  /// The id of the layer the cluster computes on.
  [[nodiscard]] std::size_t layerId() const;
  [[nodiscard]] HaloType getLayerType() const;

  protected:
  /// A cluster on the layer `layer`, named `kind` in its description, e.g. "compute" or "dr";
  /// `faces` says whether it works on faces, for the queue of its stream.
  ComputeCluster(double maxTimeStepSize,
                 long timeStepRate,
                 Executor executor,
                 initializer::Layer* layer,
                 std::string kind,
                 bool faces,
                 unsigned int profilingId,
                 LoopStatistics* loopStatistics,
                 DeviceTimeline* deviceTimeline,
                 ActorStateStatistics* actorStateStatistics);

  ActorEvent recordActionEvent() override;
  void waitForEvent(const ActorEvent& event) override;
  void timeSet(double time) override;

  initializer::Layer* layer_;
  seissol::parallel::runtime::StreamRuntime streamRuntime_;
  ClusterClock clock_;
  LoopStatistics* loopStatistics_;
  DeviceTimer deviceLoops_;
  ActorStateStatistics* actorStateStatistics_;
  unsigned int profilingId_;
  bool runTimeOutputs_{false};
#ifdef ACL_DEVICE
  device::DeviceInstance& device_ = device::DeviceInstance::instance();
#endif

  private:
  std::string kind_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_CORE_TIMESTEPPING_COMPUTE_COMPUTECLUSTER_H_

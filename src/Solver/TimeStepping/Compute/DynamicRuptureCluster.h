// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_DYNAMICRUPTURECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_DYNAMICRUPTURECLUSTER_H_

#include "Common/Executor.h"
#include "Common/Real.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "DynamicRupture/Output/PickpointSampling.h"
#include "Initializer/Typedefs.h"
#include "Kernels/DynamicRupture.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/DeviceTimer.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Actor/StepParams.h"
#include "Solver/TimeStepping/Compute/ClusterClock.h"
#include "Solver/TimeStepping/Compute/FaceCluster.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::solver {

/**
 * What the time manager sees of a dynamic rupture cluster, whatever the configuration of its faces.
 */
class DynamicRuptureClusterInterface : public FaceCluster {
  public:
  using FaceCluster::FaceCluster;

  virtual void setFaultOutputManager(dr::output::OutputManager* outputManager) = 0;
};

/**
 * The dynamic rupture faces of one layer, all of the configuration `Cfg`: the space-time
 * interpolation of the adjacent cells, the friction law, and the on-fault receiver output.
 */
template <typename Cfg>
class DynamicRuptureCluster : public DynamicRuptureClusterInterface {
  public:
  using real = Real<Cfg>;

  DynamicRuptureCluster(double maxTimeStepSize,
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
                        ActorStateStatistics* actorStateStatistics);

  ~DynamicRuptureCluster() override = default;

  ActResult act() override;

  void finalize() override;
  void finishPhase() override;

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;

  /// the copy of the clock the device writes back after each correction (see ClusterClock)
  [[nodiscard]] const double* deviceTime() const override {
    return isDeviceOn() ? clock_.host() : nullptr;
  }

  void setFaultOutputManager(dr::output::OutputManager* outputManager) override {
    faultOutputManager_ = outputManager;
  }

  [[nodiscard]] std::size_t layerId() const;

  [[nodiscard]] HaloType getLayerType() const;

  [[nodiscard]] std::string description() const override;

  protected:
  void interact(const StepParams& params) override;
  ActorEvent recordActionEvent() override;
  void waitForEvent(void* event) override;
  void timeSet(double time) override;
  StepWork prepare(ActorAction action) override;

  public:
  [[nodiscard]] bool outputsAhead(long steps) const override;
  [[nodiscard]] bool hostWork() const override;
  [[nodiscard]] void* outputStream() override;
  void setRunTimeOutputs(bool runTimeOutputs) override;

  protected:
  private:
  void computeDynamicRupture(const StepParams& params);
  void computeDynamicRuptureDevice(const StepParams& params);
  /// decides at which times of the current step the fault receivers record
  void planPickpointOutput(const StepParams& params);
  /// records the fault receivers at the times `times` of the step; `lastStep` tells whether it is
  /// the last step before the synchronization point
  void writePickpointOutput(const StepParams& params,
                            const std::vector<double>& times,
                            bool lastStep);

  /// decides at which times of the step that starts at `stepTime` the fault receivers record, and
  /// records them; `lastStep` tells whether it is the last step before the synchronization point
  void recordPickpointsNow(double stepTime, bool lastStep);

  //! the fault receivers decide about their samples when the work of a step runs
  bool runTimeOutputs_{false};
  //! the start of the step the fault receivers record in, for decisions at run time, and its
  //! address on the device, which writes it
  double pickpointStepStartHost_{0};
  double* pickpointStepStart_{&pickpointStepStartHost_};
  double* pickpointStepStartOnDevice_{nullptr};

  //! the times of the current step at which the fault receivers record
  std::vector<double> pickpointTimes_;
  PerformanceEstimate computeFlops();

  seissol::SeisSol& seissolInstance_;
  seissol::parallel::runtime::StreamRuntime streamRuntime_;

  //! the start of the current step, for the work on the stream
  ClusterClock clock_;
  kernels::DynamicRupture<Cfg> dynamicRuptureKernel_;

  CompoundGlobalData<Cfg> globalData_;
#ifdef ACL_DEVICE
  device::DeviceInstance& device_ = device::DeviceInstance::instance();
#endif

  DynamicRupture::Layer* layerData_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolver_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolverDevice_;
  dr::output::OutputManager* faultOutputManager_;

  //! when the on-fault receivers take their samples
  dr::output::PickpointSampling pickpointSampling_;

  PerformanceEstimate estimate_{};
  std::size_t perfHandle_{};

  LoopStatistics* loopStatistics_;
  //! times the loops enqueued on the device
  DeviceTimer deviceLoops_;
  ActorStateStatistics* actorStateStatistics_;
  unsigned regionComputeDynamicRupture_{};
  //! the name of the loop in a trace (see monitoring::Trace)
  std::uint32_t traceDynamicRupture_{};

  //! id used to identify this cluster when profiling
  unsigned int profilingId_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_COMPUTE_DYNAMICRUPTURECLUSTER_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTER_H_

#include "Core/Common/Executor.h"
#include "Core/Common/Real.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/DeviceTimer.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ClusterClock.h"
#include "Core/TimeStepping/Compute/FaceCluster.h"
#include "Features/DynamicRupture/Output/PickpointOutput.h"
#include "Features/DynamicRupture/Output/PickpointSampling.h"
#include "Features/DynamicRupture/Solver/DynamicRupture.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

/// An output of receivers on the fault, and when its receivers take their samples.
struct PickpointPlan {
  dr::output::PickpointOutput* output{nullptr};
  dr::output::PickpointSampling sampling;
};

/**
 * The dynamic rupture faces of one layer, all of the configuration `Cfg`: the space-time
 * interpolation of the adjacent cells, the friction law, and the on-fault receiver output.
 */
template <typename Cfg>
class DynamicRuptureCluster : public FaceCluster {
  public:
  using real = Real<Cfg>;

  DynamicRuptureCluster(double maxTimeStepSize,
                        long timeStepRate,
                        Executor executor,
                        unsigned int profilingId,
                        const std::vector<PickpointPlan>& pickpoints,
                        CompoundGlobalData<Cfg> globalData,
                        DynamicRupture::Layer* layerData,
                        const dr::friction_law::FrictionSolverFactory& frictionSolverFactory,
                        const dr::friction_law::FrictionSolverFactory& frictionSolverFactoryDevice,
                        monitoring::FlopCounter& flopCounter,
                        LoopStatistics* loopStatistics,
                        DeviceTimeline* deviceTimeline,
                        ActorStateStatistics* actorStateStatistics);

  ~DynamicRuptureCluster() override = default;

  void finalize() override;

  protected:
  void interact(const StepParams& params) override;
  StepWork prepare(ActorAction action) override;

  public:
  [[nodiscard]] bool outputsAhead(long steps) const override;
  [[nodiscard]] bool hostWork() const override;
  /// The work of one step of the faces as the FLOP counter counts it: that of their space-time
  /// interpolation (the friction law counts none).
  [[nodiscard]] PerformanceEstimate interactionEstimate() const { return estimate_; }
  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() override;
  void setRunTimeOutputs(bool runTimeOutputs) override;

  protected:
  private:
  void computeDynamicRupture(const StepParams& params);
  void computeDynamicRuptureDevice(const StepParams& params);
  /// decides at which times of the current step the receivers of each output on the fault record
  void planPickpointOutput(const StepParams& params);
  /// the times planPickpointOutput() decided, one list per output
  [[nodiscard]] std::vector<std::vector<double>> plannedPickpointTimes() const;
  /// records the receivers on the fault at the times `times` of the step, a list per output;
  /// `lastStep` tells whether it is the last step before the synchronization point
  void writePickpointOutput(const StepParams& params,
                            const std::vector<std::vector<double>>& times,
                            bool lastStep);

  /// decides at which times of the step that starts at `stepTime` the receivers of the output
  /// `slot` record, and records them; `lastStep` tells whether it is the last step before the
  /// synchronization point
  void recordPickpointsNow(std::size_t slot, double stepTime, bool lastStep);

  //! an output of receivers on the fault that the cluster samples
  struct PickpointSlot {
    dr::output::PickpointOutput* output{nullptr};
    //! when its receivers take their samples
    dr::output::PickpointSampling sampling;
    //! the times of the current step at which its receivers record
    std::vector<double> times;
    //! the start of the step its receivers record in, for decisions at run time (in pinned memory
    //! with a device, otherwise in stepStartHost), and its address on the device, which writes it
    double stepStartHost{0};
    double* stepStart{nullptr};
    double* stepStartOnDevice{nullptr};
  };
  //! the outputs, in the order of the parameters; they get their receivers when the output gets
  //! set up, after the clusters, so the cluster asks for those of its layer when it samples
  std::vector<PickpointSlot> pickpoints_;

  PerformanceEstimate computeFlops();

  monitoring::FlopCounter& flopCounter_;
  kernels::DynamicRupture<Cfg> dynamicRuptureKernel_;

  CompoundGlobalData<Cfg> globalData_;

  DynamicRupture::Layer* layerData_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolver_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolverDevice_;

  PerformanceEstimate estimate_{};
  std::size_t perfHandle_{};

  unsigned regionComputeDynamicRupture_{};
  //! the name of the loop in a trace (see monitoring::Trace)
  std::uint32_t traceDynamicRupture_{};
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTER_H_

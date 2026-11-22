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
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/ActorStateStatistics.h"
#include "Core/Monitoring/DeviceTimer.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Actor/ActorState.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/ClusterClock.h"
#include "Core/TimeStepping/Compute/FaceCluster.h"
#include "Features/DynamicRupture/Output/OutputManager.h"
#include "Features/DynamicRupture/Output/PickpointSampling.h"
#include "Features/DynamicRupture/Solver/DynamicRupture.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {} // namespace seissol

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
                        RunContext& run,
                        LoopStatistics* loopStatistics,
                        DeviceTimeline* deviceTimeline,
                        ActorStateStatistics* actorStateStatistics);

  ~DynamicRuptureCluster() override = default;

  void finalize() override;

  void setFaultOutputManager(dr::output::OutputManager* outputManager) override {
    faultOutputManager_ = outputManager;
  }

  protected:
  void interact(const StepParams& params) override;
  StepWork prepare(ActorAction action) override;

  public:
  [[nodiscard]] bool outputsAhead(long steps) const override;
  [[nodiscard]] bool hostWork() const override;
  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() override;
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

  //! the start of the step the fault receivers record in, for decisions at run time, and its
  //! address on the device, which writes it
  double pickpointStepStartHost_{0};
  double* pickpointStepStart_{&pickpointStepStartHost_};
  double* pickpointStepStartOnDevice_{nullptr};

  //! the times of the current step at which the fault receivers record
  std::vector<double> pickpointTimes_;
  PerformanceEstimate computeFlops();

  RunContext& run_;
  kernels::DynamicRupture<Cfg> dynamicRuptureKernel_;

  CompoundGlobalData<Cfg> globalData_;

  DynamicRupture::Layer* layerData_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolver_;
  std::unique_ptr<dr::friction_law::FrictionSolverImpl<Cfg>> frictionSolverDevice_;
  dr::output::OutputManager* faultOutputManager_;

  //! when the on-fault receivers take their samples
  dr::output::PickpointSampling pickpointSampling_;

  PerformanceEstimate estimate_{};
  std::size_t perfHandle_{};

  unsigned regionComputeDynamicRupture_{};
  //! the name of the loop in a trace (see monitoring::Trace)
  std::uint32_t traceDynamicRupture_{};
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SOLVER_DYNAMICRUPTURECLUSTER_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_DYNAMICRUPTURECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_DYNAMICRUPTURECLUSTER_H_

#include "Common/Executor.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/DynamicRupture.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/LoopStatistics.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/ActorState.h"
#include "Solver/TimeStepping/FaceCluster.h"
#include "Solver/TimeStepping/StepParams.h"

#include <cstddef>
#include <memory>
#include <string>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::time_stepping {

/**
 * The dynamic rupture faces of one layer: the space-time interpolation of the adjacent cells, the
 * friction law, and the on-fault receiver output.
 */
class DynamicRuptureCluster : public FaceCluster {
  public:
  DynamicRuptureCluster(double maxTimeStepSize,
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
                        ActorStateStatistics* actorStateStatistics);

  ~DynamicRuptureCluster() override = default;

  ActResult act() override;

  void finalize() override;

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;

  void setFaultOutputManager(dr::output::OutputManager* outputManager) {
    faultOutputManager_ = outputManager;
  }

  [[nodiscard]] std::size_t layerId() const;

  [[nodiscard]] HaloType getLayerType() const;

  [[nodiscard]] std::string description() const override;

  protected:
  void interact(const StepParams& params) override;
  void* recordActionEvent() override;
  void waitForEvent(void* event) override;

  private:
  void computeDynamicRupture(const StepParams& params);
  void computeDynamicRuptureDevice(const StepParams& params);
  void writePickpointOutput(const StepParams& params);
  PerformanceEstimate computeFlops();

  seissol::SeisSol& seissolInstance_;
  seissol::parallel::runtime::StreamRuntime streamRuntime_;
  kernels::DynamicRupture dynamicRuptureKernel_;

  GlobalData* globalDataOnHost_{nullptr};
  GlobalData* globalDataOnDevice_{nullptr};
#ifdef ACL_DEVICE
  device::DeviceInstance& device_ = device::DeviceInstance::getInstance();
#endif

  DynamicRupture::Layer* layerData_;
  std::unique_ptr<dr::friction_law::FrictionSolver> frictionSolver_;
  std::unique_ptr<dr::friction_law::FrictionSolver> frictionSolverDevice_;
  dr::output::OutputManager* faultOutputManager_;

  //! time interval of the on-fault receiver output
  double outputTimestep_;

  PerformanceEstimate estimate_{};
  std::size_t perfHandle_{};

  LoopStatistics* loopStatistics_;
  ActorStateStatistics* actorStateStatistics_;
  unsigned regionComputeDynamicRupture_{};

  //! id used to identify this cluster when profiling
  unsigned int profilingId_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_DYNAMICRUPTURECLUSTER_H_

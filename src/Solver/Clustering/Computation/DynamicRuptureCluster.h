// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FACECLUSTER_H_
#define FACECLUSTER_H_

#include "Kernels/DynamicRupture.h"
#include "Solver/Clustering/AbstractTimeCluster.h"
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <DynamicRupture/Output/OutputManager.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/Tree/Layer.h>
#include <Monitoring/ActorStateStatistics.h>
#include <Parallel/Helper.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/Clustering/ActorState.h>
#include <memory>

namespace seissol::solver::clustering::computation {
class DynamicRuptureCluster : public FaceCluster {
  public:
  DynamicRuptureCluster(double maxTimeStepSize,
                        long timeStepRate,
                        unsigned profilingId,
                        LayerType layerType,
                        seissol::initializer::Layer* layer,
                        seissol::initializer::DynamicRupture* descr,
                        CompoundGlobalData globalData,
                        seissol::dr::friction_law::FrictionSolver* frictionSolver,
                        seissol::dr::friction_law::FrictionSolver* frictionSolverDevice,
                        seissol::dr::output::OutputManager* faultOutputManager,
                        seissol::SeisSol& seissolInstance,
                        LoopStatistics* loopStatistics,
                        ActorStateStatistics* actorStateStatistics,
                        const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor);

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream) override;
  void computeFlops();

  void setFaultOutputManager(seissol::dr::output::OutputManager* faultOutputManager) {
    // TODO(David): remove
    this->faultOutputManager = faultOutputManager;
  }
  std::string description() const override {
    if (getLayerType() == Interior) {
      return "interior-face";
    } else {
      return "copy-face";
    }
  }
  ~DynamicRuptureCluster() override = default;

  LayerType getLayerType() const override { return layerType; }

  protected:
  seissol::SeisSol& seissolInstance;

  LayerType layerType;

  void computeDynamicRupture();
  void handleDynamicRupture();
  void computeDynamicRuptureDevice();

  void computeDynamicRuptureFlops(long long& nonZeroFlops, long long& hardwareFlops);

  void handleAdvancedComputeTimeMessage(ComputeStep step,
                                        const NeighborCluster& neighborCluster) override {}
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override {}

  seissol::initializer::Layer* layer;
  seissol::initializer::DynamicRupture* descr;
  seissol::dr::friction_law::FrictionSolver* frictionSolver;
  seissol::dr::friction_law::FrictionSolver* frictionSolverDevice;
  std::unique_ptr<seissol::dr::friction_law::FrictionSolver> frictionSolverDeviceStorage{nullptr};
  seissol::dr::output::OutputManager* faultOutputManager;

  LoopStatistics* loopStatistics;
  ActorStateStatistics* actorStateStatistics;
  unsigned regionComputeDynamicRupture;

  seissol::kernels::DynamicRupture dynamicRuptureKernel;

  GlobalData* globalDataOnHost{nullptr};
  GlobalData* globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  bool computePickpoint{false};

  ClusterTimes ct;

  long long flopsNonZero{0};
  long long flopsHardware{0};

  unsigned profilingId;

  void start() override {}

  void runCompute(ComputeStep step) override;

  private:
  void initFrictionSolverDevice();
};
} // namespace seissol::solver::clustering::computation

#endif

// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_WAVECLUSTER_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_WAVECLUSTER_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Real.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/DeviceTimeline.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Actor/ActorStateStatistics.h"
#include "Core/TimeStepping/Actor/StepParams.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Compute/CellStage.h"
#include "WavePropagation/Output/Receiver.h"
#include "WavePropagation/Physics/InitialField.h"
#include "WavePropagation/Solver/ConfigBoundary.h"
#include "WavePropagation/Solver/Solver.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

namespace seissol::solver {

/**
 * The cells of one layer, all of the configuration `Cfg`, whose waves the discontinuous Galerkin
 * scheme propagates: a prediction integrates them locally in time and space, a correction adds
 * the fluxes from their neighbors, and runs the cell stages of the features on each cell after
 * it. The cells also sample the receivers that lie in them.
 **/
template <typename Cfg>
class WaveCluster : public CellCluster {
  public:
  using real = Real<Cfg>;

  /// What a wave cluster takes from the run besides its cells.
  struct Setup {
    RunContext& run;
    //! the global data of all configurations, for the faces between configurations on the device
    GlobalDataStore& globalData;
    //! the configurations of the run
    std::vector<ConfigId> configs;
    //! the initial conditions of the cells of the configuration `Cfg`
    const std::vector<std::unique_ptr<physics::InitialField>>* initialConditions{nullptr};
    double gravitationalAcceleration{0};
  };

  WaveCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              double maxTimeStepSize,
              long timeStepRate,
              Executor executor,
              bool printProgress,
              LTS::Layer* clusterData,
              const Setup& setup,
              LoopStatistics* loopStatistics,
              DeviceTimeline* deviceTimeline,
              ActorStateStatistics* actorStateStatistics);

  ~WaveCluster() override = default;

  void addCellStage(std::unique_ptr<CellStageBase> stage) override;

  void addReceiverCluster(kernels::ReceiverCluster* receiverCluster) override;

  void finalize() override;

  [[nodiscard]] bool outputsAhead(long steps) const override;
  [[nodiscard]] std::vector<parallel::runtime::StreamRuntime*> runtimes() override;
  void setTime(double time) override;
  void finishPhase() override;

  protected:
  std::function<void()> hostPrediction(const StepParams& params) override;
  void devicePrediction(const StepParams& params) override;
  void hostCorrection(const StepParams& params) override;
  void deviceCorrection(const StepParams& params) override;
  bool planPrediction(const StepParams& params) override;
  void planCorrection(const StepParams& params) override;
  void synchronizeShared(initializer::AllocationPlace place) override;

  private:
  /**
   * The local integration of the cells: their time integration, their volume integral and their
   * local boundary integral. Their DOFs lack the contributions of their neighbors afterwards.
   **/
  void computeLocalIntegration(const StepParams& params);
  void computeLocalIntegrationDevice(const StepParams& params);

  /**
   * The contributions of the neighboring cells to the boundary integrals, after which the DOFs are
   * at the next time step (without those of dynamic rupture faces); the cell stages run on each
   * cell right after its contributions.
   **/
  void computeNeighboringIntegration(const StepParams& params);
  void computeNeighboringIntegrationDevice(const StepParams& params);

  /**
   * Samples the receivers of the cluster, with the samplings `samplings` planned for the step (one
   * per output, in the order of the outputs).
   **/
  void writeReceivers(const StepParams& params,
                      const std::vector<kernels::ReceiverCluster::Sampling>& samplings);

  /// the samplings of the receivers planned for the current step, in the order of the outputs
  [[nodiscard]] std::vector<kernels::ReceiverCluster::Sampling> plannedSamplings() const;

  void computeLocalIntegrationFlops();
  void computeNeighborIntegrationFlops();
  void computeFlops();

  enum class ComputePart : std::size_t { Local = 0, Neighbor, DRNeighbor, NumComputeParts };
  void incrementPerformanceMetrics(ComputePart part);

  RunContext& run_;

  kernels::Spacetime<Cfg> spacetimeKernel_;
  kernels::Time<Cfg> timeKernel_;
  kernels::Local<Cfg> localKernel_;
  kernels::Neighbor<Cfg> neighborKernel_;

  //! the neighbors of the cells that compute in other configurations
  kernels::ConfigBoundary<Cfg> configBoundary_;

  CompoundGlobalData<Cfg> globalData_;
  double gravitationalAcceleration_;

  //! the work features add to the corrections over the cells, in the order it runs in; at most
  //! MaxCellStages of them, whose counts the loop over the cells sums up without allocating
  static constexpr std::size_t MaxCellStages = 4;
  std::vector<std::unique_ptr<CellStage<Cfg>>> cellStages_;

  std::array<PerformanceEstimate, static_cast<std::size_t>(ComputePart::NumComputeParts)>
      estimate_{};
  std::array<std::size_t, static_cast<std::size_t>(ComputePart::NumComputeParts)> perfHandle_{};

  unsigned regionComputeLocalIntegration_{};
  unsigned regionComputeNeighboringIntegration_{};
  //! the names of the loops in a trace (see monitoring::Trace)
  std::uint32_t traceLocalIntegration_{};
  std::uint32_t traceNeighboringIntegration_{};

  //! the receivers of an output the cluster samples, the time of their next sample, and their
  //! samples in the current step, as prepare() plans them
  struct ReceiverSlot {
    kernels::ReceiverCluster* cluster{nullptr};
    double nextTime{0};
    kernels::ReceiverCluster::Sampling sampling{};
  };
  std::vector<ReceiverSlot> receivers_;
  //! the time the receivers start from (see setTime())
  double receiverStartTime_{0};
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_WAVECLUSTER_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Executor.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Monitoring/FlopCounter.h"
#include "Core/Monitoring/LoopStatistics.h"
#include "Core/Monitoring/Metric.h"
#include "Core/TimeStepping/Compute/ComputeCluster.h"
#include "Features/DynamicRupture/Solver/DynamicRuptureCluster.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Features/Plasticity/Solver/Plasticity.h"
#include "Features/Plasticity/Solver/PlasticityStage.h"
#include "Simulation/Bench/LoopBench.h"
#include "Simulation/Bench/Proxy.h"
#include "Simulation/Bench/SyntheticData.h"
#include "WavePropagation/Solver/WaveCluster.h"

#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef LIKWID_PERFMON
#include "Core/Common/Instrumentation.h"
#endif

namespace seissol::proxy {

namespace {

/**
 * The loops of the clusters of a simulation in the configuration `Cfg`, on SyntheticData: a
 * WaveCluster of the cells, with the plasticity as a stage of its correction if the cells have
 * one, and a DynamicRuptureCluster of the fault faces with their friction law. They run step by
 * step, outside of the time stepping.
 */
template <typename Cfg>
class ClusterBench : public LoopBench {
  public:
  explicit ClusterBench(const BenchSetup& setup);
  ~ClusterBench() override;

  ClusterBench(const ClusterBench&) = delete;
  ClusterBench(ClusterBench&&) = delete;
  auto operator=(const ClusterBench&) = delete;
  auto operator=(ClusterBench&&) = delete;

  void run(Kernel kernel) override;
  void wait() override;
  [[nodiscard]] auto estimate(Kernel kernel) const -> PerformanceEstimate override;

  private:
  /// The cluster of the faults; throws without faults.
  [[nodiscard]] solver::DynamicRuptureCluster<Cfg>& faultCluster() const;

  /// `cluster`, once the work of the cluster that ran before, if another one, is done: each reads
  /// what the other one writes.
  template <typename ClusterT>
  ClusterT& use(ClusterT& cluster) {
    if (last_ != nullptr && last_ != &cluster) {
      waitFor(*last_);
    }
    last_ = &cluster;
    return cluster;
  }

  static void waitFor(solver::ComputeCluster& cluster);
  void waitForAll();

  /// Runs `work`, marked as the region `name` for LIKWID on the host.
  template <typename F>
  void marked(const char* name, F&& work);

  Executor executor_;
  SyntheticData<Cfg> data_;
  monitoring::FlopCounter flopCounter_;
  LoopStatistics loopStatistics_;
  //! the work of the plasticity in a correction
  PerformanceEstimate plasticityEstimate_{};

  std::unique_ptr<solver::WaveCluster<Cfg>> waveCluster_;
  std::unique_ptr<solver::DynamicRuptureCluster<Cfg>> faultCluster_;
  solver::ComputeCluster* last_{nullptr};
};

template <typename Cfg>
ClusterBench<Cfg>::ClusterBench(const BenchSetup& setup) : executor_(setup.executor), data_(setup) {
  // the regions of the loops of the clusters, as the time manager adds them
  loopStatistics_.addRegion("computeLocalIntegration");
  loopStatistics_.addRegion("computeNeighboringIntegration");
  loopStatistics_.addRegion("computeDynamicRupture");

  LTS::Layer& layer = data_.cells();
  const typename solver::WaveCluster<Cfg>::Setup waveSetup{.flopCounter = flopCounter_,
                                                           .globalData = data_.globalData(),
                                                           .configs = {configIdOf<Cfg>()},
                                                           .initialConditions = nullptr,
                                                           .gravitationalAcceleration =
                                                               GravitationalAcceleration};
  waveCluster_ = std::make_unique<solver::WaveCluster<Cfg>>(/*clusterId=*/0,
                                                            /*globalClusterId=*/0,
                                                            layer.id(),
                                                            Timestep,
                                                            /*timeStepRate=*/1,
                                                            executor_,
                                                            /*printProgress=*/false,
                                                            &layer,
                                                            waveSetup,
                                                            &loopStatistics_,
                                                            nullptr,
                                                            nullptr);

  if (data_.plasticity() != CellPlasticity::None) {
    waveCluster_->addCellStage(std::make_unique<kernels::PlasticityStage<Cfg>>(
        flopCounter_, RelaxationTime, data_.globalData().template get<Cfg>(), layer));
    // every cell is checked, and yields with the data that make it yield
    const auto [check, yield] = kernels::Plasticity<Cfg>::metrics();
    plasticityEstimate_ = check * layer.size();
    if (data_.plasticity() == CellPlasticity::Yield) {
      plasticityEstimate_ += yield * layer.size();
    }
  }

  if (data_.hasFaces()) {
    DynamicRupture::Layer& faces = data_.faces();
    faultCluster_ =
        std::make_unique<solver::DynamicRuptureCluster<Cfg>>(Timestep,
                                                             /*timeStepRate=*/1,
                                                             executor_,
                                                             faces.id(),
                                                             std::vector<solver::PickpointPlan>{},
                                                             data_.globalData().template get<Cfg>(),
                                                             &faces,
                                                             data_.frictionLaw(),
                                                             data_.frictionLawDevice(),
                                                             flopCounter_,
                                                             &loopStatistics_,
                                                             nullptr,
                                                             nullptr);
  }
}

template <typename Cfg>
ClusterBench<Cfg>::~ClusterBench() {
  // as at the end of a run, the clusters give their streams and events back
  waitForAll();
  if (faultCluster_ != nullptr) {
    faultCluster_->finalize();
  }
  waveCluster_->finalize();
}

template <typename Cfg>
solver::DynamicRuptureCluster<Cfg>& ClusterBench<Cfg>::faultCluster() const {
  if (faultCluster_ == nullptr) {
    throw std::invalid_argument("The loops of the faces of dynamic rupture need the faces.");
  }
  return *faultCluster_;
}

template <typename Cfg>
void ClusterBench<Cfg>::waitFor(solver::ComputeCluster& cluster) {
  for (auto* runtime : cluster.runtimes()) {
    runtime->wait();
  }
}

template <typename Cfg>
void ClusterBench<Cfg>::waitForAll() {
  waitFor(*waveCluster_);
  if (faultCluster_ != nullptr) {
    waitFor(*faultCluster_);
  }
}

template <typename Cfg>
template <typename F>
void ClusterBench<Cfg>::marked([[maybe_unused]] const char* name, F&& work) {
#ifdef LIKWID_PERFMON
  const bool host = executor_ == Executor::Host;
  if (host) {
#pragma omp parallel
    {
      LIKWID_MARKER_START(name);
    }
  }
#endif
  std::invoke(std::forward<F>(work));
#ifdef LIKWID_PERFMON
  if (host) {
#pragma omp parallel
    {
      LIKWID_MARKER_STOP(name);
    }
  }
#endif
}

template <typename Cfg>
void ClusterBench<Cfg>::run(Kernel kernel) {
  switch (kernel) {
  case Kernel::Local:
    marked("local", [&] { use(*waveCluster_).predictOnce(Step); });
    break;
  case Kernel::Neighbor:
  case Kernel::NeighborDR:
    marked("neighboring", [&] { use(*waveCluster_).correctOnce(Step); });
    break;
  case Kernel::GodunovDR:
    use(faultCluster()).interactOnce(Step);
    break;
  case Kernel::All:
    run(Kernel::Local);
    run(Kernel::Neighbor);
    break;
  case Kernel::AllDR:
    run(Kernel::Local);
    run(Kernel::GodunovDR);
    run(Kernel::NeighborDR);
    break;
  default:
    throw std::invalid_argument("Only the kernel mode runs the kernel " +
                                Aux::kernel2str({kernel}) + ".");
  }
}

template <typename Cfg>
void ClusterBench<Cfg>::wait() {
  waitForAll();
}

template <typename Cfg>
auto ClusterBench<Cfg>::estimate(Kernel kernel) const -> PerformanceEstimate {
  switch (kernel) {
  case Kernel::Local:
    return waveCluster_->predictionEstimate();
  case Kernel::Neighbor:
  case Kernel::NeighborDR:
    return waveCluster_->correctionEstimate() + plasticityEstimate_;
  case Kernel::GodunovDR:
    return faultCluster().interactionEstimate();
  case Kernel::All:
    return estimate(Kernel::Local) + estimate(Kernel::Neighbor);
  case Kernel::AllDR:
    return estimate(Kernel::Local) + estimate(Kernel::GodunovDR) + estimate(Kernel::NeighborDR);
  default:
    throw std::invalid_argument("Only the kernel mode runs the kernel " +
                                Aux::kernel2str({kernel}) + ".");
  }
}

} // namespace

auto makeClusterBench(const BenchSetup& setup) -> std::unique_ptr<LoopBench> {
  return dispatchConfig(setup.config, [&](auto cfg) -> std::unique_ptr<LoopBench> {
    using Cfg = decltype(cfg);
    return std::make_unique<ClusterBench<Cfg>>(setup);
  });
}

} // namespace seissol::proxy

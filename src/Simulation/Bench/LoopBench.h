// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_BENCH_LOOPBENCH_H_
#define SEISSOL_SRC_SIMULATION_BENCH_LOOPBENCH_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Monitoring/Metric.h"
#include "Simulation/Bench/Proxy.h"

#include <cstddef>
#include <memory>

namespace seissol::proxy {

/// What the loops of the proxy run on (see SyntheticData).
struct BenchSetup {
  ConfigId config{};
  std::size_t cells{};
  //! with faces of dynamic rupture: all faces of the cells lie on the fault
  bool dynamicRupture{false};
  CellPlasticity plasticity{CellPlasticity::None};
  //! (with dynamic rupture)
  FrictionLaw frictionLaw{};
  Executor executor{Executor::Host};
};

/**
 * The loops of a time step of a simulation on synthetic data (see SyntheticData), one step at a
 * time: kernel by kernel (KernelBench), or as the clusters of a simulation run them (ClusterBench).
 */
class LoopBench {
  public:
  LoopBench() = default;
  virtual ~LoopBench() = default;

  LoopBench(const LoopBench&) = delete;
  LoopBench(LoopBench&&) = delete;
  auto operator=(const LoopBench&) = delete;
  auto operator=(LoopBench&&) = delete;

  /// Runs the loops of `kernel` once; on the device, they may still run afterwards (see wait()).
  virtual void run(Kernel kernel) = 0;

  /// Waits until the loops run so far are done.
  virtual void wait() = 0;

  /// The work of one run of the loops of `kernel`, as a simulation counts it.
  [[nodiscard]] virtual auto estimate(Kernel kernel) const -> PerformanceEstimate = 0;
};

/// The loops of the mode `mode` on the data `setup`.
auto makeLoopBench(Mode mode, const BenchSetup& setup) -> std::unique_ptr<LoopBench>;

/// The kernels of the kernel mode, and the clusters of the cluster mode, on the data `setup`.
auto makeKernelBench(const BenchSetup& setup) -> std::unique_ptr<LoopBench>;
auto makeClusterBench(const BenchSetup& setup) -> std::unique_ptr<LoopBench>;

} // namespace seissol::proxy

#endif // SEISSOL_SRC_SIMULATION_BENCH_LOOPBENCH_H_

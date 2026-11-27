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

/**
 * The loops of a time step of a simulation, on synthetic data: an interior layer of cells of one
 * configuration, with random neighbors and arbitrary values. Either all faces of the cells are
 * regular, or all of them lie on the fault, with four times as many fault faces as cells, whose
 * adjacent cells are random as well.
 *
 * The loops are those of the clusters of a simulation (see WaveCluster and DynamicRuptureCluster),
 * run step by step outside of the time stepping: every cell with the same time step, of size 1,
 * and no friction on the fault. The cells read the time integrals of their neighbors from the
 * buffers of these (GTS), and no feature works on them.
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

/// The loops on `cells` cells of the configuration `config`, with faces on the fault if
/// `dynamicRupture`, running on `executor`.
auto makeLoopBench(ConfigId config, std::size_t cells, bool dynamicRupture, Executor executor)
    -> std::unique_ptr<LoopBench>;

} // namespace seissol::proxy

#endif // SEISSOL_SRC_SIMULATION_BENCH_LOOPBENCH_H_

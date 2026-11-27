// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_BENCH_PROXY_H_
#define SEISSOL_SRC_SIMULATION_BENCH_PROXY_H_

#include "Alignment.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Monitoring/Metric.h"

#include <cstdint>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace seissol::proxy {

/**
 * The loops of a time step that the proxy measures, on the data of LoopBench; one after the other
 * in each step, as a list of them names them. In both modes:
 *
 *   local        the time and the local integration of the cells
 *   neigh        the contributions of their neighbors (neigh_dr: with all faces on the fault)
 *   godunov_dr   the work of the fault faces
 *   all          local, neigh: a time step without dynamic rupture
 *   all_dr       local, godunov_dr, neigh_dr: a time step with dynamic rupture
 *
 * In the kernel mode only, the parts of these:
 *
 *   ader         the time integration of the cells (the first part of local)
 *   localwoader  their local integration (the second part of local)
 *   plast        the plastic correction of cells that stay within the yield surface
 *   plast_yield  the plastic correction of cells that yield
 */
enum class Kernel {
  All = 0,
  Local,
  Neighbor,
  Ader,
  LocalWOAder,
  NeighborDR,
  GodunovDR,
  AllDR,
  Plasticity,
  PlasticityYield
};

/**
 * How the proxy runs the loops. The kernel mode runs each kernel over all cells or faces on its
 * own, as the proxy did before it ran the clusters; godunov_dr is the space-time interpolation of
 * the faces alone there. The cluster mode runs the loops of the clusters of a simulation
 * (WaveCluster, DynamicRuptureCluster), with all they do around the kernels: godunov_dr is the
 * work of a cluster of fault faces there, the friction law included, and the correction of the
 * cells runs the plasticity with them if they have one.
 */
enum class Mode { Kernel, Cluster };

/// The plasticity of the cells in the cluster mode; in the kernel mode, the kernels plast and
/// plast_yield choose it.
enum class CellPlasticity : std::uint8_t {
  None,
  //! cells that stay within the yield surface: the correction only checks them
  Elastic,
  //! cells that yield in every step
  Yield
};

enum class OutputFormat { Plain, Json };

struct ProxyConfig {
  std::uint64_t cells{static_cast<std::uint64_t>(1e5)};
  std::uint64_t timesteps{10};
  //! the loops of a step, in the order they run in
  std::vector<Kernel> kernels;
  Mode mode{Mode::Kernel};
  //! (in the cluster mode)
  CellPlasticity plasticity{CellPlasticity::None};
  bool verbose{true};
  //! where the loops run: on the host, or (in a GPU build, by default) on the device
  Executor executor{isDeviceOn() ? Executor::Device : Executor::Host};
  //! The configuration the kernels run in.
  ConfigId config{defaultConfig()};
  //! For how many seconds the kernels run before the timed steps (besides one step in any case): a
  //! device raises its clocks only under load, within some hundred milliseconds.
  double warmup{0};
};

struct ProxyOutput {
  double time{};
  double cycles{};
  /// Which hardware counter `cycles` came from; see Simulation/Bench/Cycles.h.
  std::string_view cycleSource{"none"};
  //! "kernel" or "cluster", see Mode
  std::string_view mode{"kernel"};
  double libxsmmNumTotalGFlop{};
  double pspammNumTotalGFlop{};
  double libxsmmAndpspammNumTotalGFlop{};
  double actualNonZeroGFlop{};
  double actualHardwareGFlop{};
  double gib{};
  double kernelGib{};
  double nonZeroFlopPerCycle{};
  double hardwareFlopPerCycle{};
  double bytesPerCycle{};
  double kernelBytesPerCycle{};
  double nonZeroGFlops{};
  double hardwareGFlops{};
  double gibPerSecond{};
  double kernelGibPerSecond{};
};

struct Aux {
  static auto kernel2str(const std::vector<Kernel>& kernels) -> std::string;

  static auto str2kernel(const std::string& kernelStr) -> std::vector<Kernel>;

  static auto getAllowedKernels() -> std::vector<std::string>;

  static auto mode2str(Mode mode) -> std::string_view;

  static void writeOutput(std::ostream& stream,
                          const ProxyOutput& output,
                          const std::string& kernelStr,
                          OutputFormat format);

  static void displayOutput(const ProxyOutput& output, const std::string& kernelStr);
};

/// Whether one of the loops `kernels` works on the faces of dynamic rupture.
auto needsDynamicRupture(const std::vector<Kernel>& kernels) -> bool;

/// The plasticity the cells of `config` need.
auto cellPlasticity(const ProxyConfig& config) -> CellPlasticity;

/// Throws std::invalid_argument if the loops of `config` do not go together, e.g. a kernel of
/// the kernel mode in the cluster mode.
void checkConfig(const ProxyConfig& config);

/// Times `config.timesteps` steps of the loops `config.kernels` on `config.cells` cells.
auto runProxy(const ProxyConfig& config) -> ProxyOutput;

/// The work of one time step of the kernels of `config` on its cells, as the proxy counts it,
/// without running them.
auto estimateProxy(const ProxyConfig& config) -> PerformanceEstimate;

} // namespace seissol::proxy

#endif // SEISSOL_SRC_SIMULATION_BENCH_PROXY_H_

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

#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace seissol::proxy {

/**
 * The loops of a time step that the proxy measures, on the data of LoopBench: the prediction of
 * the cells (local), their correction (neigh, or neigh_dr with all faces of the cells on the
 * fault), the work of the fault faces (godunov_dr: the space-time interpolation, and the friction
 * law of a fault without friction), and the loops of a time step without dynamic rupture (all:
 * local, neigh) and with it (all_dr: local, godunov_dr, neigh_dr).
 */
enum class Kernel { All = 0, Local, Neighbor, NeighborDR, GodunovDR, AllDR };

enum class OutputFormat { Plain, Json };

struct ProxyConfig {
  unsigned cells{static_cast<unsigned>(1e5)};
  unsigned timesteps{10};
  //! the loops of a step, in the order they run in
  std::vector<Kernel> kernels;
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

  static void writeOutput(std::ostream& stream,
                          const ProxyOutput& output,
                          const std::string& kernelStr,
                          OutputFormat format);

  static void displayOutput(const ProxyOutput& output, const std::string& kernelStr);
};

/// Whether one of the loops `kernels` works on the faces of dynamic rupture.
auto needsDynamicRupture(const std::vector<Kernel>& kernels) -> bool;

/// Times `config.timesteps` steps of the loops `config.kernels` on `config.cells` cells.
auto runProxy(const ProxyConfig& config) -> ProxyOutput;

/// The work of one time step of the kernels of `config` on its cells, as the proxy counts it,
/// without running them.
auto estimateProxy(const ProxyConfig& config) -> PerformanceEstimate;

} // namespace seissol::proxy

#endif // SEISSOL_SRC_SIMULATION_BENCH_PROXY_H_

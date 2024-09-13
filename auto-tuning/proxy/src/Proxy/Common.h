// SPDX-FileCopyrightInfo: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_PROXY_PROXY_COMMON_DATATYPES_HPP
#define SEISSOL_PROXY_PROXY_COMMON_DATATYPES_HPP

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::proxy {

enum class Kernel { All = 0, Local, Neighbor, Ader, LocalWOAder, NeighborDR, GodunovDR, AllDR };

enum class OutputFormat { Plain, Json };

struct ProxyConfig {
  unsigned cells{static_cast<unsigned>(1e5)};
  unsigned timesteps{10};
  Kernel kernel{Kernel::All};
  bool verbose{true};
};

struct ProxyOutput {
  double time{};
  double cycles{};
  double libxsmmNumTotalGFlop{};
  double pspammNumTotalGFlop{};
  double libxsmmAndpspammNumTotalGFlop{};
  double actualNonZeroGFlop{};
  double actualHardwareGFlop{};
  double gib{};
  double nonZeroFlopPerCycle{};
  double hardwareFlopPerCycle{};
  double bytesPerCycle{};
  double nonZeroGFlops{};
  double hardwareGFlops{};
  double gibPerSecond{};
};

ProxyOutput runProxy(ProxyConfig config);

struct Aux {
  static auto kernel2str(Kernel kernel) -> std::string;

  static auto str2kernel(const std::string& kernelStr) -> Kernel;

  static auto getAllowedKernels() -> std::vector<std::string>;

  static void writeOutput(std::ostream& stream,
                          const ProxyOutput& output,
                          const std::string& kernelStr,
                          OutputFormat format);

  static void displayOutput(const ProxyOutput& output, const std::string& kernelStr);
};

} // namespace seissol::proxy

#endif // SEISSOL_PROXY_PROXY_COMMON_DATATYPES_HPP

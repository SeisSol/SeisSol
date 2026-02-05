// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PROXY_COMMON_H_
#define SEISSOL_SRC_PROXY_COMMON_H_

#include "Common/Executor.h"

#include <cstdint>
#include <string>
#include <vector>

namespace seissol::proxy {

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
  PlasticityYield,
  FL0,
  FL3,
  FL3TP,
  FL4,
  FL4TP,
  FL6,
  FL15,
  FL16,
  FL33,
  FL34,
  FL103,
  FL103TP
};

enum class OutputFormat : std::uint32_t { Plain = 0, Json, PlainTFlop };

struct ProxyConfig {
  std::uint64_t cells{static_cast<std::uint64_t>(1e6)};
  std::uint64_t timesteps{10};
  std::vector<Kernel> kernels;
  bool verbose{true};
  Executor executor;
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

} // namespace seissol::proxy

#endif // SEISSOL_SRC_PROXY_COMMON_H_

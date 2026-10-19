// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common.h"

#include "Common/Json.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utils/stringutils.h>
#include <vector>

namespace {
using namespace seissol::proxy;
const std::unordered_map<Kernel, std::string> Map{
    {Kernel::All, "all"},
    {Kernel::Local, "local"},
    {Kernel::Neighbor, "neigh"},
    {Kernel::Ader, "ader"},
    {Kernel::LocalWOAder, "localwoader"},
    {Kernel::NeighborDR, "neigh_dr"},
    {Kernel::GodunovDR, "godunov_dr"},
    {Kernel::AllDR, "all_dr"},
};

const std::unordered_map<std::string, Kernel> InvMap{{"all", Kernel::All},
                                                     {"local", Kernel::Local},
                                                     {"neigh", Kernel::Neighbor},
                                                     {"ader", Kernel::Ader},
                                                     {"localwoader", Kernel::LocalWOAder},
                                                     {"neigh_dr", Kernel::NeighborDR},
                                                     {"godunov_dr", Kernel::GodunovDR},
                                                     {"all_dr", Kernel::AllDR}};
} // namespace

namespace seissol::proxy {

auto Aux::kernel2str(const std::vector<Kernel>& kernel) -> std::string {
  std::string output;
  for (std::size_t i = 0; i < kernel.size(); ++i) {
    if (i > 0) {
      output += ",";
    }
    if (Map.find(kernel[i]) != Map.end()) {
      output += Map.at(kernel[i]);
    } else {
      throw std::runtime_error("Unknown kernel type");
    }
  }
  return output;
}

auto Aux::str2kernel(const std::string& kernelStr) -> std::vector<Kernel> {
  auto parts = utils::StringUtils::split(kernelStr, ',');
  for (auto& part : parts) {
    utils::StringUtils::trim(part);
  }

  std::vector<Kernel> kernels(parts.size());
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (InvMap.find(parts[i]) != InvMap.end()) {
      kernels[i] = InvMap.at(parts[i]);
    } else {
      throw std::runtime_error("Unknown kernel string: " + parts[i]);
    }
  }

  return kernels;
}

auto Aux::getAllowedKernels() -> std::vector<std::string> {
  std::vector<std::string> kernels{};
  kernels.reserve(InvMap.size());
  for (const auto& kernel : InvMap) {
    kernels.push_back(kernel.first);
  }
  return kernels;
}

void Aux::writeOutput(std::ostream& stream,
                      const ProxyOutput& output,
                      const std::string& kernelStr,
                      OutputFormat format) {
  if (format == OutputFormat::Plain) {
    stream << "=================================================\n";
    stream << "===            PERFORMANCE SUMMARY            ===\n";
    stream << "=================================================\n";
    stream << "seissol proxy mode                  : " << kernelStr << '\n';
    stream << "time for seissol proxy              : " << output.time << '\n';
    stream << "cycles                              : " << output.cycles << '\n';
    stream << '\n';
    stream << "GFLOP (libxsmm)                     : " << output.libxsmmNumTotalGFlop << '\n';
    stream << "GFLOP (pspamm)                      : " << output.pspammNumTotalGFlop << '\n';
    stream << "GFLOP (libxsmm + pspamm)            : " << output.libxsmmAndpspammNumTotalGFlop
           << '\n';
    stream << "GFLOP (non-zero) for seissol proxy  : " << output.actualNonZeroGFlop << '\n';
    stream << "GFLOP (hardware) for seissol proxy  : " << output.actualHardwareGFlop << '\n';
    stream << "GiB (estimate) for seissol proxy    : " << output.gib << '\n';
    stream << "GiB (kernels) for seissol proxy     : " << output.kernelGib << '\n';
    stream << '\n';
    stream << "FLOPS/cycle (non-zero)              : " << output.nonZeroFlopPerCycle << '\n';
    stream << "FLOPS/cycle (hardware)              : " << output.hardwareFlopPerCycle << '\n';
    stream << "Bytes/cycle (estimate)              : " << output.bytesPerCycle << '\n';
    stream << "Bytes/cycle (kernel)                : " << output.kernelBytesPerCycle << '\n';
    stream << '\n';
    stream << "GFLOPS (non-zero) for seissol proxy : " << output.nonZeroGFlops << '\n';
    stream << "GFLOPS (hardware) for seissol proxy : " << output.hardwareGFlops << '\n';
    stream << "GiB/s (estimate) for seissol proxy  : " << output.gibPerSecond << '\n';
    stream << "GiB/s (kernel) for seissol proxy    : " << output.kernelGibPerSecond << '\n';
    stream << "=================================================\n";
    stream << '\n';
  } else {
    // Written for machine consumption, on one line. No non-finite value may reach it: JSON has no
    // spelling for those, and the readers expect a number for every field. Runner.cpp keeps every
    // field finite; this is the second line of defence.
    const auto finite = [](double value) { return std::isfinite(value) ? value : 0.0; };
    json::Value document;
    document["name"] = kernelStr;
    document["cycle-source"] = output.cycleSource;
    document["time"] = finite(output.time);
    document["cycles"] = finite(static_cast<double>(output.cycles));
    document["gflop-libxsmm"] = finite(output.libxsmmNumTotalGFlop);
    document["gflop-pspamm"] = finite(output.pspammNumTotalGFlop);
    document["gflop-libxsmm-pspamm"] = finite(output.libxsmmAndpspammNumTotalGFlop);
    document["gflop-nz"] = finite(output.actualNonZeroGFlop);
    document["gflop-hw"] = finite(output.actualHardwareGFlop);
    document["gib"] = finite(output.gib);
    document["gib-kernel"] = finite(output.kernelGib);
    document["gflopcycle-nz"] = finite(output.nonZeroFlopPerCycle);
    document["gflopcycle-hw"] = finite(output.hardwareFlopPerCycle);
    document["gibcycle"] = finite(output.bytesPerCycle);
    document["gibcycle-kernel"] = finite(output.kernelBytesPerCycle);
    document["gflops-nz"] = finite(output.nonZeroGFlops);
    document["gflops-hw"] = finite(output.hardwareGFlops);
    document["gibs"] = finite(output.gibPerSecond);
    document["gibs-kernel"] = finite(output.kernelGibPerSecond);
    stream << document.dump(-1) << '\n';
  }
}

void Aux::displayOutput(const ProxyOutput& output, const std::string& kernelStr) {
  writeOutput(std::cout, output, kernelStr, OutputFormat::Plain);
}

} // namespace seissol::proxy

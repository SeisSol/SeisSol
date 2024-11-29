#ifndef SEISSOL_PROXY_PROXY_COMMON_DATATYPES_HPP
#define SEISSOL_PROXY_PROXY_COMMON_DATATYPES_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

enum class Kernel {
  All = 0,
  Local,
  Neighbor,
  Ader,
  LocalWOAder,
  NeighborDR,
  GodunovDR,
  Dynrup
};

struct ProxyConfig {
  unsigned cells{static_cast<unsigned>(1e5)};
  unsigned timesteps{10};
  unsigned phase{10};
  unsigned fault{0};
  Kernel kernel{Kernel::All};
  bool verbose{true};
};

struct ProxyOutput{
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
  static std::string kernel2str(Kernel kernel) {
    if (map.find(kernel) != map.end()) {
      return map[kernel];
    } else {
      throw std::runtime_error("unknown kernel type");
    }
  }

  static Kernel str2kernel(const std::string& kernelStr) {
    if (invMap.find(kernelStr) != invMap.end()) {
      return invMap[kernelStr];
    } else {
      throw std::runtime_error("unknown kernel string: " + kernelStr);
    }
  }

  static std::vector<std::string> getAllowedKernels() {
    std::vector<std::string> kernels{};
    kernels.reserve(invMap.size());
    for (const auto& kernel: invMap) {
      kernels.push_back(kernel.first);
    }
    return kernels;
  }

  static void displayOutput(const ProxyOutput& output, const std::string& kernelStr) {
    printf("=================================================\n");
    printf("===            PERFORMANCE SUMMARY            ===\n");
    printf("=================================================\n");
    printf("seissol proxy mode                  : %s\n",   kernelStr.c_str());
    printf("time for seissol proxy              : %f\n",   output.time);
    printf("cycles                              : %f\n\n", output.cycles);
    printf("GFLOP (libxsmm)                     : %f\n",   output.libxsmmNumTotalGFlop);
    printf("GFLOP (pspamm)                      : %f\n",   output.pspammNumTotalGFlop);
    printf("GFLOP (libxsmm + pspamm)            : %f\n",   output.libxsmmAndpspammNumTotalGFlop);
    printf("GFLOP (non-zero) for seissol proxy  : %f\n",   output.actualNonZeroGFlop);
    printf("GFLOP (hardware) for seissol proxy  : %f\n",   output.actualHardwareGFlop);
    printf("GiB (estimate) for seissol proxy    : %f\n\n", output.gib);
    printf("FLOPS/cycle (non-zero)              : %f\n",   output.nonZeroFlopPerCycle);
    printf("FLOPS/cycle (hardware)              : %f\n",   output.hardwareFlopPerCycle);
    printf("Bytes/cycle (estimate)              : %f\n\n", output.bytesPerCycle);
    printf("GFLOPS (non-zero) for seissol proxy : %f\n",   output.nonZeroGFlops);
    printf("GFLOPS (hardware) for seissol proxy : %f\n",   output.hardwareGFlops);
    printf("GiB/s (estimate) for seissol proxy  : %f\n",   output.gibPerSecond);
    printf("=================================================\n");
    printf("\n");
  }

protected:
  inline static std::unordered_map<Kernel, std::string> map{
      {Kernel::All,         "all"},
      {Kernel::Local,       "local"},
      {Kernel::Neighbor,       "neigh"},
      {Kernel::Ader,        "ader"},
      {Kernel::LocalWOAder, "localwoader"},
      {Kernel::NeighborDR,    "neigh_dr"},
      {Kernel::GodunovDR,  "godunov_dr"},
      {Kernel::Dynrup,  "dynrup"},
  };

  inline static std::unordered_map<std::string, Kernel> invMap{
      {"all", Kernel::All},
      {"local", Kernel::Local},
      {"neigh", Kernel::Neighbor},
      {"ader", Kernel::Ader},
      {"localwoader", Kernel::LocalWOAder},
      {"neigh_dr", Kernel::NeighborDR},
      {"godunov_dr", Kernel::GodunovDR},
      {"dynrup", Kernel::Dynrup},
  };
};


#endif //SEISSOL_PROXY_PROXY_COMMON_DATATYPES_HPP

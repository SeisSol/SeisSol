// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Alignment.h"
#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/Executor.h"
#include "Core/Common/Instrumentation.h"
#include "Simulation/Bench/Proxy.h"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <utils/args.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Core/Common/Typedefs.h"
#include "Core/Parallel/Helper.h"
#include "Core/Parallel/Runtime/Stream.h"

#include <Device/device.h>
#include <cstdlib>
#endif

using namespace seissol::proxy;

namespace {
int runMain(int argc, char* argv[]) {
  std::stringstream kernelHelp;
  const auto allowedKernels = Aux::getAllowedKernels();
  kernelHelp << "Kernels to benchmark. A comma-separated list of (those kernels will be run "
                "sequentially in each pass; ader, localwoader, plast and plast_yield in the kernel "
                "mode only): ";
  {
    bool comma = false;
    for (const auto& kernel : allowedKernels) {
      if (comma) {
        kernelHelp << ", ";
      }
      kernelHelp << kernel;
      comma = true;
    }
  }

  const std::vector<std::string> formatValues = {"plain", "json"};
  const std::vector<std::string> modeValues = {"kernel", "cluster"};
  const std::vector<std::string> plasticityValues = {"none", "elastic", "yield"};

  utils::Args args("The SeisSol proxy is used to benchmark the kernels used in the SeisSol "
                   "earthquake simulation software. This version of SeisSol proxy was built for "
                   "the following configurations:\n" +
                   seissol::describeBuiltConfigs());
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());
  args.addEnumOption("format", formatValues, 'f', "The output format", false);
  args.addEnumOption("mode",
                     modeValues,
                     'm',
                     "How the loops run: kernel by kernel, as the proxy ran them before it ran the "
                     "clusters (kernel, the default), or as the clusters of a simulation run them, "
                     "with all they do around the kernels (cluster)",
                     false);
  args.addEnumOption("plasticity",
                     plasticityValues,
                     'p',
                     "In the cluster mode, the plasticity of the cells: none (the default), "
                     "cells that stay within the yield surface (elastic), or cells that yield in "
                     "every step (yield)",
                     false);
  args.addOption("config",
                 'c',
                 "The configuration to run the kernels in, by its name (default: the one the "
                 "environment variable SEISSOL_CONFIGURATION names, or else the first one listed "
                 "above)",
                 utils::Args::Required,
                 false);

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  ProxyConfig config{};
  config.cells = args.getAdditionalArgument<std::size_t>("cells");
  config.timesteps = args.getAdditionalArgument<std::size_t>("timesteps");
  const auto kernelStr = args.getAdditionalArgument<std::string>("kernel");
  const auto formatValue = args.getArgument<std::int32_t>("format", 0);
  if (args.isSet("config")) {
    const auto configName = args.getArgument<std::string>("config");
    const auto configId = seissol::findConfig(configName);
    if (!configId.has_value()) {
      std::cerr << "The configuration " << configName << " is not built into this proxy." << '\n';
      return -1;
    }
    config.config = *configId;
  }

  const auto format = formatValue == 1 ? OutputFormat::Json : OutputFormat::Plain;
  config.mode = args.getArgument<std::int32_t>("mode", 0) == 1 ? Mode::Cluster : Mode::Kernel;
  config.plasticity = static_cast<CellPlasticity>(args.getArgument<std::int32_t>("plasticity", 0));

  try {
    config.kernels = Aux::str2kernel(kernelStr);
    checkConfig(config);
  } catch (std::exception& error) {
    std::cerr << error.what() << '\n';
    return -1;
  }

  config.executor = seissol::isDeviceOn() ? seissol::Executor::Device : seissol::Executor::Host;

  LIKWID_MARKER_INIT;
  // the regions of the loops over the cells on the host (see KernelBench and ClusterBench)
#pragma omp parallel
  {
    LIKWID_MARKER_REGISTER("ader");
    LIKWID_MARKER_REGISTER("localwoader");
    LIKWID_MARKER_REGISTER("local");
    LIKWID_MARKER_REGISTER("neighboring");
    LIKWID_MARKER_REGISTER("plasticity");
  }

#ifdef ACL_DEVICE
  // with SEISSOL_HSA_QUEUES (HIP only), the stream runtime of the proxy runs on the HSA queue that
  // a run of SeisSol gives all stream runtimes besides those of the clusters: the kernels then cost
  // what they cost on HSA queues
  const auto hsaQueues = seissol::hsaQueueCount();
  if (hsaQueues > 0 && seissol::Backend == seissol::DeviceBackend::Hip) {
    // (before HIP starts) the streams of HIP then only record graphs, on one hardware queue each
    // NOLINTNEXTLINE(concurrency-mt-unsafe,misc-include-cleaner)
    setenv("GPU_MAX_HW_QUEUES", "1", 0);
  }
  using DeviceType = ::device::DeviceInstance;
  auto& device = DeviceType::instance();
  device.api().setDevice(0);
  device.api().initialize();
  if (hsaQueues > 0 && device.api().enableHardwareQueues()) {
    seissol::parallel::runtime::HsaQueues queues{};
    queues.count = hsaQueues;
    seissol::parallel::runtime::useHsaQueues(queues);
  }
#endif

  const auto output = runProxy(config);

#ifdef ACL_DEVICE
  device.finalize();
#endif
  LIKWID_MARKER_CLOSE;

  Aux::writeOutput(std::cout, output, kernelStr, format);
  std::cout.flush();

  return 0;
}
} // namespace

int main(int argc, char* argv[]) {
  // (an error the proxy does not catch itself, e.g. one reading the environment, ends it with a
  // message instead of std::terminate)
  try {
    return runMain(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "The proxy stopped with an error: " << error.what() << '\n';
  }
  return 1;
}

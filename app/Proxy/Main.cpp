// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigRegistry.h"
#include "Common/Executor.h"
#include "Kernels/Common.h"
#include "Proxy/Common.h"
#include "Proxy/LikwidWrapper.h"
#include "Proxy/Runner.h"
#include "Proxy/Tools.h"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utils/args.h>
#include <vector>

#ifdef ACL_DEVICE
#include "Alignment.h"
#include "Common/Typedefs.h"
#include "Parallel/Helper.h"
#include "Parallel/Runtime/Stream.h"

#include <Device/device.h>
#include <cstdlib>
#endif

using namespace seissol::proxy;

namespace {
int runMain(int argc, char* argv[]) {
  std::stringstream kernelHelp;
  const auto allowedKernels = Aux::getAllowedKernels();
  kernelHelp << "Kernels to benchmark. A comma-separated list of (those kernels will be run "
                "sequentially in each pass): ";
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

  utils::Args args("The SeisSol proxy is used to benchmark the kernels used in the SeisSol "
                   "earthquake simulation software. This version of SeisSol proxy was built for "
                   "the following configurations:\n" +
                   seissol::describeBuiltConfigs());
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());
  args.addEnumOption("format", formatValues, 'f', "The output format", false);
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

  try {
    config.kernels = Aux::str2kernel(kernelStr);
  } catch (std::runtime_error& error) {
    std::cerr << error.what() << '\n';
    return -1;
  }

  config.executor = seissol::isDeviceOn() ? seissol::Executor::Device : seissol::Executor::Host;

  LIKWID_MARKER_INIT;

  registerMarkers();

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
  printHostname();

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

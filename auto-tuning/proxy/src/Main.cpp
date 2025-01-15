// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

#include "Proxy/Common.h"
#include "Proxy/Runner.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utils/args.h>

using namespace seissol::proxy;

auto main(int argc, char* argv[]) -> int {
  std::stringstream kernelHelp;
  auto allowedKernels = Aux::getAllowedKernels();
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

  const char* formatValues[] = {"plain", "json"};

  utils::Args args("The SeisSol proxy is used to benchmark the kernels used in the SeisSol "
                   "earthquake simulation software.");
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());
  args.addEnumOption("format", formatValues, 'f', "The output format", false);

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  ProxyConfig config{};
  config.cells = args.getAdditionalArgument<unsigned>("cells");
  config.timesteps = args.getAdditionalArgument<unsigned>("timesteps");
  const auto kernelStr = args.getAdditionalArgument<std::string>("kernel");
  const auto formatValue = args.getArgument<int>("format", 0);

  const auto format = formatValue == 1 ? OutputFormat::Json : OutputFormat::Plain;

  try {
    config.kernels = Aux::str2kernel(kernelStr);
  } catch (std::runtime_error& error) {
    std::cerr << error.what() << std::endl;
    return -1;
  }

  auto output = runProxy(config);
  Aux::writeOutput(std::cout, output, kernelStr, format);
  std::cout.flush();

  return 0;
}

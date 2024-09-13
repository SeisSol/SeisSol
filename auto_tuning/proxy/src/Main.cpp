#include "Proxy/Common.hpp"
#include <iostream>
#include <utils/args.h>

using namespace seissol::proxy;

auto main(int argc, char* argv[]) -> int {
  std::stringstream kernelHelp;
  auto allowedKernels = Aux::getAllowedKernels();
  kernelHelp << "Kernel: ";
  for (const auto& kernel : allowedKernels) {
    kernelHelp << ", " << kernel;
  }

  utils::Args args("The SeisSol proxy is used to benchmark the kernels used in the SeisSol "
                   "earthquake simulation software.");
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  ProxyConfig config{};
  config.cells = args.getAdditionalArgument<unsigned>("cells");
  config.timesteps = args.getAdditionalArgument<unsigned>("timesteps");
  auto kernelStr = args.getAdditionalArgument<std::string>("kernel");

  try {
    config.kernel = Aux::str2kernel(kernelStr);
  } catch (std::runtime_error& error) {
    std::cerr << error.what() << std::endl;
    return -1;
  }

  auto output = runProxy(config);
  Aux::writeOutput(std::cout, output, kernelStr, OutputFormat::Plain);
  return 0;
}

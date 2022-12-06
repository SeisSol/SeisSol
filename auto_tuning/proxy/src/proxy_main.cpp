#include <utils/args.h>
#include "proxy_common.hpp"
#include <iostream>


int main(int argc, char* argv[]) {
  std::stringstream kernelHelp;
  auto allowedKernels = Aux::getAllowedKernels();
  kernelHelp << "Kernel: ";
  for (const auto& kernel: allowedKernels) {
    kernelHelp << ", " << kernel;
  }

  utils::Args args;
  args.addOption("tasks", 0, "Use task based parallelization", utils::Args::No, false);
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  ProxyConfig config{};
  config.cells = args.getAdditionalArgument<unsigned>("cells");
  config.timesteps = args.getAdditionalArgument<unsigned>("timesteps");
  config.kernelConfig.parallelizationStrategy =
      args.isSet("tasks") ? ParallelizationStrategy::Taskloop : ParallelizationStrategy::ParallelFor;
  auto kernelStr = args.getAdditionalArgument<std::string>("kernel");

  try {
    config.kernel = Aux::str2kernel(kernelStr);
  }
  catch (std::runtime_error& error) {
    std::cerr << error.what() << std::endl;
    return -1;
  }

  auto output = runProxy(config);
  Aux::displayOutput(output, kernelStr);
  return 0;
}

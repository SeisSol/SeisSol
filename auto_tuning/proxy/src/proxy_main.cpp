#include <utils/args.h>
#include "proxy_common.hpp"
#include <iostream>

#include <omp.h>

void targetDartTextPointer2() {}

int main(int argc, char* argv[]) {
  // the proxy does not have MPI
#ifdef USE_TARGETDART
  initTargetDART(reinterpret_cast<void*>(&targetDartTextPointer2));
#endif

  // argument reading
  std::stringstream kernelHelp;
  auto allowedKernels = Aux::getAllowedKernels();
  kernelHelp << "Kernel: ";
  for (const auto& kernel: allowedKernels) {
    kernelHelp << ", " << kernel;
  }

  utils::Args args;
  args.addAdditionalOption("cells", "Number of cells");
  args.addAdditionalOption("timesteps", "Number of timesteps");
  args.addAdditionalOption("kernel", kernelHelp.str());
  args.addAdditionalOption("phase", "Length of synchronization phase (in timesteps)");
  args.addAdditionalOption("fault", "Fault type (only used if kernel == \"dynrup\")", false);

  if (args.parse(argc, argv) != utils::Args::Success) {
    return -1;
  }

  ProxyConfig config{};
  config.cells = args.getAdditionalArgument<unsigned>("cells");
  config.timesteps = args.getAdditionalArgument<unsigned>("timesteps");
  config.phase = args.getAdditionalArgument<unsigned>("phase");
  config.fault = args.getAdditionalArgument<unsigned>("fault", 0);
  auto kernelStr = args.getAdditionalArgument<std::string>("kernel");

  try {
    config.kernel = Aux::str2kernel(kernelStr);
  }
  catch (std::runtime_error& error) {
    std::cerr << error.what() << std::endl;
    return -1;
  }

  // proxy execution
  auto output = runProxy(config);

  // postprocessing
  Aux::displayOutput(output, kernelStr);

  // again, no MPI
#ifdef USE_TARGETDART
  finalizeTargetDART();
#endif
  return 0;
}

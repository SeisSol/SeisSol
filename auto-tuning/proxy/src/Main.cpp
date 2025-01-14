// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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

  const std::vector<std::string> formatValues = {"plain", "json"};

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

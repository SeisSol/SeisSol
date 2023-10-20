// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "Args.hpp"

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <ostream>

namespace seissol {

Args ArgParser::parseArgs(int argc, char** argv) {
  Args a = {};
  a.memType = MemoryType::Host;
  a.synchronizationTime = 0.1;
  a.showHelp = false;

  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') {
      auto const fail = [&]() {
        throw std::runtime_error("Error: unrecognized argument " + std::string(argv[i]));
      };
      if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
        a.showHelp = true;
      } else if (i + 1 < argc) {
        if (std::strcmp(argv[i], "-m") == 0 || std::strcmp(argv[i], "--mem-type") == 0) {
          ++i;
          if (std::strcmp(argv[i], "host") == 0) {
            a.memType = MemoryType::Host;
          }
#ifdef HAVE_SYCL
          else if (std::strcmp(argv[i], "shared") == 0) {
            a.memType = MemoryType::Shared;
          } else if (std::strcmp(argv[i], "device") == 0) {
            a.memType = MemoryType::Device;
          }
#endif
          else {
            fail();
          }

        } else if (std::strcmp(argv[i], "-s") == 0 || std::strcmp(argv[i], "--sync-time") == 0) {
          a.synchronizationTime = atof(argv[++i]);
        } else {
          fail();
        }
      } else {
        fail();
      }
    } else {
      a.prefix = std::string(argv[i]);
    }
  }

  return a;
}

void ArgParser::showHelp(std::ostream& os) {
  os << "usage: SeisSol_comm_proxy prefix" << std::endl
     << R"HELP(
Reproduces comm pattern of SeisSol.
Needs a JSON dump of the LTS structure that is obtained by setting
DumpLocalTimeSteppingStructure=1
in the parameters file.

positional arguments:
    prefix              Path to JSON dump (WITHOUT .X.json sufix)

optional arguments:
    -h, --help          Show help and quit
    -s, --sync-time     Synchronization time [default: 0.1]
    -m, --mem-type      Memory allocation type [host, shared, device, default: host]
)HELP";
}

} // namespace seissol

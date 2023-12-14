// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ARGS_20231019_HPP
#define ARGS_20231019_HPP

#include "MemoryType.hpp"

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace seissol {

struct Args {
  MemoryType memType;
  std::string prefix;
  double synchronizationTime;
  bool showHelp;
};

class ArgParser {
  public:
  static Args parseArgs(int argc, char** argv);
  static void showHelp(std::ostream& os);
};

} // namespace seissol

#endif // ARGS_20231019_HPP

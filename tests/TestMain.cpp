// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include "Parallel/MPI.h"

#include <string_view>

// NOLINTNEXTLINE
extern long long libxsmm_num_total_flops;
// NOLINTNEXTLINE
extern long long pspamm_num_total_flops;

namespace {
// Whether doctest is asked to list the test cases rather than to run them.
bool listsTestCases(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg(argv[i]);
    if (arg.find("list-test-cases") != std::string_view::npos || arg == "-ltc" ||
        arg == "-dt-ltc" || arg.rfind("-ltc=", 0) == 0 || arg.rfind("-dt-ltc=", 0) == 0) {
      return true;
    }
  }
  return false;
}
} // namespace

int main(int argc, char** argv) {
  // make sure these two variables are always included into the tests
  // (sometimes not the case for single-module test binaries)
  libxsmm_num_total_flops = 0;
  pspamm_num_total_flops = 0;

  seissol::Mpi::mpi.init(argc, argv);
  doctest::Context context;

  context.applyCommandLine(argc, argv);

  // CTest learns the tests from this list. Under an MPI launcher every rank would print it, which
  // registers each test once per rank -- and a launcher that interleaves the output of the ranks
  // garbles it. So only the first rank lists them.
  if (seissol::Mpi::mpi.rank() != 0 && listsTestCases(argc, argv)) {
    context.setOption("out", "/dev/null");
  }

  const int returnValue = context.run();

  seissol::Mpi::finalize();

  return returnValue;
}

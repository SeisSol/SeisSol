// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include "Parallel/MPI.h"

// NOLINTNEXTLINE
extern long long libxsmm_num_total_flops;
// NOLINTNEXTLINE
extern long long pspamm_num_total_flops;

int main(int argc, char** argv) {
  // make sure these two variables are always included into the tests
  // (sometimes not the case for single-module test binaries)
  libxsmm_num_total_flops = 0;
  pspamm_num_total_flops = 0;

  seissol::Mpi::mpi.init(argc, argv);
  doctest::Context context;

  context.applyCommandLine(argc, argv);

  const int returnValue = context.run();

  seissol::Mpi::finalize();

  return returnValue;
}

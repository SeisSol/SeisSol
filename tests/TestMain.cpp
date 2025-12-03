// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest.h>

#include "Parallel/MPI.h"

int main(int argc, char** argv) {
  seissol::Mpi::mpi.init(argc, argv);
  doctest::Context context;

  context.applyCommandLine(argc, argv);

  const int returnValue = context.run();

  seissol::Mpi::finalize();

  return returnValue;
}

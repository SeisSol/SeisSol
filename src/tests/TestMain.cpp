// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#define DOCTEST_CONFIG_IMPLEMENT
#include "Parallel/MPI.h"
#include "doctest.h"

int main(int argc, char** argv) {
  seissol::MPI::mpi.init(argc, argv);
  doctest::Context context;

  context.applyCommandLine(argc, argv);

  const int returnValue = context.run();

  seissol::MPI::finalize();

  return returnValue;
}

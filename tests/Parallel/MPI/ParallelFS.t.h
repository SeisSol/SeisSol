// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/MPI.h"
#include "Parallel/ParallelFS.h"

#include <filesystem>

namespace seissol::unit_test {

TEST_CASE("createDirectoryGlobally creates a directory that every rank sees" *
          doctest::test_suite("mpi")) {
  const auto& mpi = seissol::Mpi::mpi;
  const auto base = std::filesystem::path("parallelfs-test-directory");
  const auto nested = base / "with" / "parents";
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(base);
  }

  seissol::createDirectoryGlobally(nested);
  CHECK(std::filesystem::is_directory(nested));

  // a directory that exists already is no error
  seissol::createDirectoryGlobally(nested);
  CHECK(std::filesystem::is_directory(nested));

  mpi.barrier();
  if (mpi.rank() == 0) {
    std::filesystem::remove_all(base);
  }
}

} // namespace seissol::unit_test

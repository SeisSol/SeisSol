// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Parallel/MPI.h"
#include "Parallel/MPIDeviceSupport.h"

#include <optional>

namespace seissol::unit_test {

TEST_CASE("The transfer mode follows what the MPI library tells about device memory") {
  using Mode = Mpi::DataTransferMode;
  using parallel::chooseDataTransferMode;
  const std::optional<bool> notTold;

  SUBCASE("a mode asked for stays") {
    for (const auto support : {notTold, std::optional<bool>{true}, std::optional<bool>{false}}) {
      for (const bool unified : {false, true}) {
        CHECK(chooseDataTransferMode(Mode::Direct, support, unified) == Mode::Direct);
        CHECK(chooseDataTransferMode(Mode::CopyInCopyOutHost, support, unified) ==
              Mode::CopyInCopyOutHost);
      }
    }
  }

  SUBCASE("an MPI library that can, or does not tell, gets the device memory") {
    CHECK(chooseDataTransferMode(std::nullopt, true, false) == Mode::Direct);
    CHECK(chooseDataTransferMode(std::nullopt, notTold, false) == Mode::Direct);
  }

  SUBCASE("one that cannot gets host memory, unless the buffers are unified memory") {
    CHECK(chooseDataTransferMode(std::nullopt, false, false) == Mode::CopyInCopyOutHost);
    CHECK(chooseDataTransferMode(std::nullopt, false, true) == Mode::Direct);
  }
}

} // namespace seissol::unit_test

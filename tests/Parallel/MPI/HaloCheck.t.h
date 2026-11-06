// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// The check of the halo exchange before the time stepping, through MPI between the ranks: every
// rank sends two regions around a ring, one to each neighbor.

#include <doctest.h>

#include "Common/Real.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCheck.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/Mpi/MpiHaloTransport.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::unit_test {

namespace halochecktest {

/**
 * Checks the exchange of two regions around the ring: region 0 (in double precision) goes to the
 * right neighbor and region 1 (in single precision, large, so that MPI usually does not send it
 * eagerly) to the left one. With `oneMore`, each ghost region of region 1 expects one value more
 * than its sender sends.
 */
inline solver::HaloCheckResult checkRing(bool persistent, bool oneMore) {
  const int rank = seissol::Mpi::mpi.rank();
  const int size = seissol::Mpi::mpi.size();
  constexpr std::size_t SizeF64 = 1000;
  constexpr std::size_t SizeF32 = 40000;
  const int right = (rank + 1) % size;
  const int left = (rank + size - 1) % size;

  std::vector<double> copyF64(SizeF64, -1.0);
  std::vector<double> ghostF64(SizeF64, -1.0);
  std::vector<float> copyF32(SizeF32, -1.0F);
  std::vector<float> ghostF32(SizeF32 + (oneMore ? 1 : 0), -1.0F);
  solver::RemoteClusterPair regions;
  regions.copy.emplace_back(copyF64.data(), SizeF64, RealType::F64, right, 0);
  regions.copy.emplace_back(copyF32.data(), SizeF32, RealType::F32, left, 1);
  regions.ghost.emplace_back(ghostF64.data(), SizeF64, RealType::F64, left, 0);
  regions.ghost.emplace_back(ghostF32.data(), ghostF32.size(), RealType::F32, right, 1);

  solver::MpiHaloTransport transport(regions, persistent);
  const auto result = solver::checkHaloExchange({{&transport, &regions, "the ring"}},
                                                rank,
                                                solver::HaloMemory::Host,
                                                std::chrono::seconds(60),
                                                true);
  transport.finalize();

  const auto zero = [](auto value) { return value == 0; };
  CHECK(std::all_of(copyF64.begin(), copyF64.end(), zero));
  CHECK(std::all_of(ghostF64.begin(), ghostF64.end(), zero));
  CHECK(std::all_of(copyF32.begin(), copyF32.end(), zero));
  CHECK(std::all_of(ghostF32.begin(), ghostF32.end(), zero));
  return result;
}

} // namespace halochecktest

TEST_CASE("The check of the halo exchange passes MPI between the ranks" *
          doctest::test_suite("mpi")) {
  for (const bool persistent : {true, false}) {
    CAPTURE(persistent);
    const auto result = halochecktest::checkRing(persistent, false);
    CHECK(result.regions == 2);
    CHECK(result.bytes == 1000 * sizeof(double) + 40000 * sizeof(float));
    CHECK(result.wrongRegions == 0);
  }
}

TEST_CASE("The check of the halo exchange finds a message that is too short" *
          doctest::test_suite("mpi")) {
  const auto result = halochecktest::checkRing(true, true);
  CHECK(result.wrongRegions == 1);
  // the value beyond the message stays zero
  CHECK(result.firstWrong.find("with tag 1 has 1 of 40001 values (f32) wrong, the first at "
                               "position 40000: 0 instead of") != std::string::npos);
}

} // namespace seissol::unit_test

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Monitoring/FlopCounter.h"

namespace seissol::unit_test {
using seissol::monitoring::FlopCounter;

TEST_CASE("FlopCounter increment operations" * doctest::test_suite("monitoring")) {
  // Note: We cannot test init() or print*() here as they require MPI.
  // We test that the increment functions don't crash and that the object
  // is constructible. The actual accumulation is internal (private members),
  // so we verify the functions are callable without error.

  FlopCounter counter;

  SUBCASE("Increment local flops") {
    counter.incrementNonZeroFlopsLocal(100);
    counter.incrementHardwareFlopsLocal(150);
    counter.incrementNonZeroFlopsLocal(200);
    counter.incrementHardwareFlopsLocal(250);
    CHECK(true);
  }

  SUBCASE("Increment neighbor flops") {
    counter.incrementNonZeroFlopsNeighbor(1000);
    counter.incrementHardwareFlopsNeighbor(1500);
    CHECK(true);
  }

  SUBCASE("Increment other flops") {
    counter.incrementNonZeroFlopsOther(500);
    counter.incrementHardwareFlopsOther(700);
    CHECK(true);
  }

  SUBCASE("Increment dynamic rupture flops") {
    counter.incrementNonZeroFlopsDynamicRupture(2000);
    counter.incrementHardwareFlopsDynamicRupture(3000);
    CHECK(true);
  }

  SUBCASE("Increment plasticity flops") {
    counter.incrementNonZeroFlopsPlasticity(400);
    counter.incrementHardwareFlopsPlasticity(600);
    CHECK(true);
  }

  SUBCASE("Zero increment is valid") {
    counter.incrementNonZeroFlopsLocal(0);
    counter.incrementHardwareFlopsLocal(0);
    CHECK(true);
  }

  SUBCASE("Large values") {
    counter.incrementHardwareFlopsLocal(1'000'000'000LL);
    counter.incrementHardwareFlopsLocal(1'000'000'000LL);
    counter.incrementHardwareFlopsLocal(1'000'000'000LL);
    CHECK(true);
  }
}

} // namespace seissol::unit_test

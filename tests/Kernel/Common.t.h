// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Kernels/Common.h"

namespace seissol::unit_test {
using namespace seissol::kernels;

// ---------------------------------------------------------------------------
// getNumberOfBasisFunctions
// ---------------------------------------------------------------------------

TEST_CASE("getNumberOfBasisFunctions" * doctest::test_suite("kernel")) {
  // Formula: O*(O+1)*(O+2)/6
  CHECK(getNumberOfBasisFunctions(1) == 1);   // 1*2*3/6
  CHECK(getNumberOfBasisFunctions(2) == 4);   // 2*3*4/6
  CHECK(getNumberOfBasisFunctions(3) == 10);  // 3*4*5/6
  CHECK(getNumberOfBasisFunctions(4) == 20);  // 4*5*6/6
  CHECK(getNumberOfBasisFunctions(5) == 35);  // 5*6*7/6
  CHECK(getNumberOfBasisFunctions(6) == 56);  // 6*7*8/6
  CHECK(getNumberOfBasisFunctions(7) == 84);  // 7*8*9/6
  CHECK(getNumberOfBasisFunctions(8) == 120); // 8*9*10/6
}

// ---------------------------------------------------------------------------
// getNumberOfAlignedReals
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("getNumberOfAlignedReals" * doctest::test_suite("kernel"),
                   RealT,
                   float,
                   double) {
  SUBCASE("Already aligned") {
    // If numberOfReals * sizeof(RealT) is already a multiple of alignment,
    // no padding needed.
    const unsigned alignment = sizeof(RealT);
    unsigned n = 10;
    CHECK(getNumberOfAlignedReals<RealT>(n, alignment) == n);
  }

  SUBCASE("Padding is applied") {
    // For a larger alignment, the result should be >= input
    unsigned n = 7;
    unsigned result = getNumberOfAlignedReals<RealT>(n, Vectorsize);
    CHECK(result >= n);
    // result * sizeof(RealT) should be a multiple of the alignment
    CHECK((result * sizeof(RealT)) % Vectorsize == 0);
  }

  SUBCASE("One real") {
    unsigned result = getNumberOfAlignedReals<RealT>(1, Vectorsize);
    CHECK(result >= 1);
    CHECK((result * sizeof(RealT)) % Vectorsize == 0);
  }

  SUBCASE("Zero reals") {
    unsigned result = getNumberOfAlignedReals<RealT>(0, Vectorsize);
    CHECK(result == 0);
  }
}

// ---------------------------------------------------------------------------
// getNumberOfAlignedBasisFunctions
// ---------------------------------------------------------------------------

TEST_CASE_TEMPLATE("getNumberOfAlignedBasisFunctions" * doctest::test_suite("kernel"),
                   RealT,
                   float,
                   double) {
  SUBCASE("At least as many as unaligned") {
    for (unsigned order = 1; order <= 8; ++order) {
      auto aligned = getNumberOfAlignedBasisFunctions<RealT>(order, Vectorsize);
      auto unaligned = getNumberOfBasisFunctions(order);
      CHECK(aligned >= unaligned);
    }
  }

  SUBCASE("Alignment property holds") {
    for (unsigned order = 1; order <= 8; ++order) {
      auto aligned = getNumberOfAlignedBasisFunctions<RealT>(order, Vectorsize);
      CHECK((aligned * sizeof(RealT)) % Vectorsize == 0);
    }
  }
}

// ---------------------------------------------------------------------------
// isDeviceOn
// ---------------------------------------------------------------------------

TEST_CASE("isDeviceOn reflects build config" * doctest::test_suite("kernel")) {
  // This just verifies the function is callable; the result depends on build config
#ifdef ACL_DEVICE
  CHECK(isDeviceOn() == true);
#else
  CHECK(isDeviceOn() == false);
#endif
}

} // namespace seissol::unit_test

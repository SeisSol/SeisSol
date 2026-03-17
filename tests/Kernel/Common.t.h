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

TEST_CASE("getNumberOfBasisFunctions default uses ConvergenceOrder" *
          doctest::test_suite("kernel")) {
  // Default argument is ConvergenceOrder (from Config)
  auto bf = getNumberOfBasisFunctions();
  auto bfExplicit = getNumberOfBasisFunctions(ConvergenceOrder);
  CHECK(bf == bfExplicit);
  CHECK(bf > 0);
}

// ---------------------------------------------------------------------------
// getNumberOfAlignedReals
// ---------------------------------------------------------------------------

TEST_CASE("getNumberOfAlignedReals" * doctest::test_suite("kernel")) {
  SUBCASE("Already aligned") {
    // If numberOfReals * sizeof(real) is already a multiple of alignment,
    // no padding needed.
    unsigned alignment = sizeof(real);
    unsigned n = 10;
    CHECK(getNumberOfAlignedReals(n, alignment) == n);
  }

  SUBCASE("Padding is applied") {
    // For a larger alignment, the result should be >= input
    unsigned n = 7;
    unsigned result = getNumberOfAlignedReals(n);
    CHECK(result >= n);
    // result * sizeof(real) should be a multiple of the alignment
    CHECK((result * sizeof(real)) % Vectorsize == 0);
  }

  SUBCASE("One real") {
    unsigned result = getNumberOfAlignedReals(1);
    CHECK(result >= 1);
    CHECK((result * sizeof(real)) % Vectorsize == 0);
  }

  SUBCASE("Zero reals") {
    unsigned result = getNumberOfAlignedReals(0);
    CHECK(result == 0);
  }
}

// ---------------------------------------------------------------------------
// getNumberOfAlignedBasisFunctions
// ---------------------------------------------------------------------------

TEST_CASE("getNumberOfAlignedBasisFunctions" * doctest::test_suite("kernel")) {
  SUBCASE("At least as many as unaligned") {
    for (unsigned order = 1; order <= 8; ++order) {
      auto aligned = getNumberOfAlignedBasisFunctions(order);
      auto unaligned = getNumberOfBasisFunctions(order);
      CHECK(aligned >= unaligned);
    }
  }

  SUBCASE("Alignment property holds") {
    for (unsigned order = 1; order <= 8; ++order) {
      auto aligned = getNumberOfAlignedBasisFunctions(order);
      CHECK((aligned * sizeof(real)) % Vectorsize == 0);
    }
  }

  SUBCASE("Default argument matches ConvergenceOrder") {
    auto a = getNumberOfAlignedBasisFunctions();
    auto b = getNumberOfAlignedBasisFunctions(ConvergenceOrder);
    CHECK(a == b);
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

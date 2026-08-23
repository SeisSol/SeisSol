// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Typedefs.h"

#include <cmath>

namespace seissol::unit_test {
using namespace seissol::dr;

// ---------------------------------------------------------------------------
// ImpedancesAndEta: physical consistency
// ---------------------------------------------------------------------------

TEST_CASE("ImpedancesAndEta default zero" * doctest::test_suite("dynamicrupture")) {
  ImpedancesAndEta imp{};
  CHECK(imp.zp == doctest::Approx(0.0));
  CHECK(imp.zs == doctest::Approx(0.0));
  CHECK(imp.zpNeig == doctest::Approx(0.0));
  CHECK(imp.zsNeig == doctest::Approx(0.0));
  CHECK(imp.etaP == doctest::Approx(0.0));
  CHECK(imp.etaS == doctest::Approx(0.0));
}

TEST_CASE("ImpedancesAndEta physical setup" * doctest::test_suite("dynamicrupture")) {
  // Typical setup from the existing FrictionSolverCommon test
  ImpedancesAndEta imp;
  imp.zp = 10.0;
  imp.zs = 20.0;
  imp.zpNeig = 15.0;
  imp.zsNeig = 25.0;
  imp.etaP = imp.zp * imp.zpNeig / (imp.zp + imp.zpNeig);
  imp.etaS = imp.zs * imp.zsNeig / (imp.zs + imp.zsNeig);
  imp.invZp = 1.0 / imp.zp;
  imp.invZs = 1.0 / imp.zs;
  imp.invZpNeig = 1.0 / imp.zpNeig;
  imp.invZsNeig = 1.0 / imp.zsNeig;

  SUBCASE("Eta is harmonic mean") {
    // etaP = zp * zpNeig / (zp + zpNeig) = 10*15/25 = 6
    CHECK(imp.etaP == doctest::Approx(6.0));
    // etaS = zs * zsNeig / (zs + zsNeig) = 20*25/45 = 500/45 ≈ 11.111
    CHECK(imp.etaS == doctest::Approx(500.0 / 45.0));
  }

  SUBCASE("Inverses are correct") {
    CHECK(imp.invZp == doctest::Approx(0.1));
    CHECK(imp.invZs == doctest::Approx(0.05));
    CHECK(imp.invZpNeig == doctest::Approx(1.0 / 15.0));
    CHECK(imp.invZsNeig == doctest::Approx(0.04));
  }

  SUBCASE("Eta between the two impedances") {
    // Harmonic mean is always <= arithmetic mean
    CHECK(imp.etaP <= (imp.zp + imp.zpNeig) / 2.0);
    CHECK(imp.etaS <= (imp.zs + imp.zsNeig) / 2.0);
    CHECK(imp.etaP > 0.0);
    CHECK(imp.etaS > 0.0);
  }
}

TEST_CASE("ImpedancesAndEta equal impedances" * doctest::test_suite("dynamicrupture")) {
  ImpedancesAndEta imp;
  imp.zp = 10.0;
  imp.zpNeig = 10.0;
  imp.etaP = imp.zp * imp.zpNeig / (imp.zp + imp.zpNeig);

  // Harmonic mean of equal values = value / 2
  CHECK(imp.etaP == doctest::Approx(5.0));
}

// ---------------------------------------------------------------------------
// ReceiverPoint: default initialization
// ---------------------------------------------------------------------------

TEST_CASE("ReceiverPoint defaults" * doctest::test_suite("dynamicrupture")) {
  ReceiverPoint rp;
  CHECK(rp.faultFaceIndex == -1);
  CHECK(rp.localFaceSideId == -1);
  CHECK(rp.elementIndex == -1);
  CHECK(rp.globalReceiverIndex == -1);
  CHECK(rp.isInside == false);
  CHECK(rp.nearestGpIndex == -1);
  CHECK(rp.faultTag == -1);
  CHECK(rp.simIndex == 0);
}

// ---------------------------------------------------------------------------
// FaultDirections: default zero
// ---------------------------------------------------------------------------

TEST_CASE("FaultDirections defaults" * doctest::test_suite("dynamicrupture")) {
  FaultDirections fd;
  for (int i = 0; i < 3; ++i) {
    CHECK(fd.faceNormal[i] == doctest::Approx(0.0));
    CHECK(fd.tangent1[i] == doctest::Approx(0.0));
    CHECK(fd.tangent2[i] == doctest::Approx(0.0));
    CHECK(fd.strike[i] == doctest::Approx(0.0));
    CHECK(fd.dip[i] == doctest::Approx(0.0));
  }
}

} // namespace seissol::unit_test

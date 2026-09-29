// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "DynamicRupture/Output/Geometry.h"
#include "DynamicRupture/Typedefs.h"

#include <cmath>
#include <cstddef>

namespace seissol::unit_test {
using namespace seissol::dr;

// ---------------------------------------------------------------------------
// ImpedancesAndEta: physical consistency
// ---------------------------------------------------------------------------

TEST_CASE("ImpedancesAndEta default zero" * doctest::test_suite("dynamicrupture")) {
  ImpedancesAndEta imp{};
  for (std::size_t point = 0; point < ImpedancePoints; ++point) {
    CHECK(imp.zp(point) == doctest::Approx(0.0));
    CHECK(imp.zs(point) == doctest::Approx(0.0));
    CHECK(imp.zpNeig(point) == doctest::Approx(0.0));
    CHECK(imp.zsNeig(point) == doctest::Approx(0.0));
    CHECK(imp.etaP(point) == doctest::Approx(0.0));
    CHECK(imp.etaS(point) == doctest::Approx(0.0));
  }
}

TEST_CASE("ImpedancesAndEta physical setup" * doctest::test_suite("dynamicrupture")) {
  // Typical setup from the existing FrictionSolverCommon test
  ImpedancesAndEta imp;
  imp.zp.fill(10.0);
  imp.zs.fill(20.0);
  imp.zpNeig.fill(15.0);
  imp.zsNeig.fill(25.0);
  imp.etaP.fill(10.0 * 15.0 / (10.0 + 15.0));
  imp.etaS.fill(20.0 * 25.0 / (20.0 + 25.0));
  imp.invZp.fill(1.0 / 10.0);
  imp.invZs.fill(1.0 / 20.0);
  imp.invZpNeig.fill(1.0 / 15.0);
  imp.invZsNeig.fill(1.0 / 25.0);

  SUBCASE("Eta is harmonic mean") {
    // etaP = zp * zpNeig / (zp + zpNeig) = 10*15/25 = 6
    CHECK(imp.etaP(0) == doctest::Approx(6.0));
    // etaS = zs * zsNeig / (zs + zsNeig) = 20*25/45 = 500/45 ≈ 11.111
    CHECK(imp.etaS(0) == doctest::Approx(500.0 / 45.0));
  }

  SUBCASE("Inverses are correct") {
    CHECK(imp.invZp(0) == doctest::Approx(0.1));
    CHECK(imp.invZs(0) == doctest::Approx(0.05));
    CHECK(imp.invZpNeig(0) == doctest::Approx(1.0 / 15.0));
    CHECK(imp.invZsNeig(0) == doctest::Approx(0.04));
  }

  SUBCASE("Eta between the two impedances") {
    // Harmonic mean is always <= arithmetic mean
    CHECK(imp.etaP(0) <= (imp.zp(0) + imp.zpNeig(0)) / 2.0);
    CHECK(imp.etaS(0) <= (imp.zs(0) + imp.zsNeig(0)) / 2.0);
    CHECK(imp.etaP(0) > 0.0);
    CHECK(imp.etaS(0) > 0.0);
  }

  SUBCASE("One value reaches every point") {
    // Where the material does not vary the scalar is stored once, and every
    // point of the face has to see it.
    for (std::size_t point = 0; point < ImpedancePoints; ++point) {
      CHECK(imp.zp(point) == doctest::Approx(10.0));
      CHECK(imp.etaS(point) == doctest::Approx(500.0 / 45.0));
    }
  }
}

TEST_CASE("ImpedancesAndEta equal impedances" * doctest::test_suite("dynamicrupture")) {
  ImpedancesAndEta imp;
  imp.zp.fill(10.0);
  imp.zpNeig.fill(10.0);
  imp.etaP.fill(imp.zp(0) * imp.zpNeig(0) / (imp.zp(0) + imp.zpNeig(0)));

  // Harmonic mean of equal values = value / 2
  CHECK(imp.etaP(0) == doctest::Approx(5.0));
}

TEST_CASE("Per-point impedances" * doctest::test_suite("dynamicrupture")) {
  // The pointwise storage has to keep the points apart, and the constant one
  // has to answer for every point from the single value it keeps. Both shapes
  // are checked here, whichever the build itself uses.
  SUBCASE("a value per point") {
    ImpedancesAndEtaOf<true> imp{};
    for (std::size_t point = 0; point < misc::NumPaddedPoints; ++point) {
      imp.zp.set(point, 1.0 + point);
    }
    for (std::size_t point = 0; point < misc::NumPaddedPoints; ++point) {
      CHECK(imp.zp(point) == doctest::Approx(1.0 + point));
    }
  }

  SUBCASE("one value for the face") {
    ImpedancesAndEtaOf<false> imp{};
    imp.zp.fill(7.0);
    for (std::size_t point = 0; point < misc::NumPaddedPoints; ++point) {
      CHECK(imp.zp(point) == doctest::Approx(7.0));
    }
    // the last write wins, since there is only one value
    imp.zp.set(misc::NumPaddedPoints - 1, 9.0);
    CHECK(imp.zp(0) == doctest::Approx(9.0));
  }

  SUBCASE("a matrix per point") {
    PointMatrix<true, 4> matrix{};
    for (std::size_t point = 0; point < misc::NumPaddedPoints; ++point) {
      matrix.at(point)[2] = static_cast<real>(point);
    }
    for (std::size_t point = 0; point < misc::NumPaddedPoints; ++point) {
      CHECK(matrix.at(point)[2] == doctest::Approx(point));
      CHECK(matrix.at(point)[0] == doctest::Approx(0.0));
    }
  }
}

// ---------------------------------------------------------------------------
// Receiver: default initialization
// ---------------------------------------------------------------------------

TEST_CASE("Receiver defaults" * doctest::test_suite("dynamicrupture")) {
  Receiver rp;
  CHECK_FALSE(rp.faultFaceIndex.hasValue());
  CHECK_FALSE(rp.localFaceSideId.hasValue());
  CHECK_FALSE(rp.localNeighborFaceSideId.hasValue());
  CHECK_FALSE(rp.elementIndex.hasValue());
  CHECK_FALSE(rp.elementGlobalIndex.hasValue());
  CHECK_FALSE(rp.elementNeighborGlobalIndex.hasValue());
  CHECK_FALSE(rp.globalReceiverIndex.hasValue());
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

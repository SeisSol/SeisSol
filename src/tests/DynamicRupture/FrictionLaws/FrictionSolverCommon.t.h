// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <numeric>

#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "tests/TestHelper.h"

#ifndef USE_ACOUSTIC

namespace seissol::unit_test {

using namespace seissol;
using namespace seissol::dr;

TEST_CASE("Friction Solver Common") {
  FaultStresses<Executor::Host> faultStresses{};
  TractionResults<Executor::Host> tractionResults{};
  ImpedancesAndEta impAndEta;
  alignas(Alignment) real qInterpolatedPlus[ConvergenceOrder][tensor::QInterpolated::size()] = {{}};
  alignas(Alignment)
      real qInterpolatedMinus[ConvergenceOrder][tensor::QInterpolated::size()] = {{}};
  alignas(Alignment) real imposedStatePlus[tensor::QInterpolated::size()] = {};
  alignas(Alignment) real imposedStateMinus[tensor::QInterpolated::size()] = {};
  double timeWeights[ConvergenceOrder];
  std::iota(std::begin(timeWeights), std::end(timeWeights), 1);
  constexpr real Epsilon = 1e4 * std::numeric_limits<real>::epsilon();

  using QInterpolatedShapeT = real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
  auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));
  using ImposedStateShapeT = real(*)[misc::NumPaddedPoints];
  auto* iSPlus = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* iSMinus = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  impAndEta.zp = 10.0;
  impAndEta.zs = 20.0;
  impAndEta.zpNeig = 15.0;
  impAndEta.zsNeig = 25.0;
  impAndEta.etaP = impAndEta.zp * impAndEta.zpNeig / (impAndEta.zp + impAndEta.zpNeig);
  impAndEta.etaS = impAndEta.zs * impAndEta.zsNeig / (impAndEta.zs + impAndEta.zsNeig);
  impAndEta.invZp = 1.0 / impAndEta.zp;
  impAndEta.invZs = 1.0 / impAndEta.zs;
  impAndEta.invZpNeig = 1.0 / impAndEta.zpNeig;
  impAndEta.invZsNeig = 1.0 / impAndEta.zsNeig;

  ImpedanceMatrices impMats;
  auto etaView = init::eta::view::create(impMats.eta);
  etaView(0, 0) = impAndEta.etaP;
  etaView(1, 1) = impAndEta.etaS;
  etaView(2, 2) = impAndEta.etaS;
  auto impedanceView = init::Zplus::view::create(impMats.impedance);
  impedanceView(0, 0) = impAndEta.invZp;
  impedanceView(1, 1) = impAndEta.invZs;
  impedanceView(2, 2) = impAndEta.invZs;
  auto impedanceNeigView = init::Zminus::view::create(impMats.impedanceNeig);
  impedanceNeigView(0, 0) = impAndEta.invZpNeig;
  impedanceNeigView(1, 1) = impAndEta.invZsNeig;
  impedanceNeigView(2, 2) = impAndEta.invZsNeig;

  auto qP = [](size_t o, size_t q, size_t p) { return static_cast<real>(o + q + p); };
  auto qM = [](size_t o, size_t q, size_t p) { return static_cast<real>(2 * (o + q + p)); };
  auto t1 = [](size_t o, size_t p) { return static_cast<real>(o + p); };
  auto t2 = [](size_t o, size_t p) { return static_cast<real>(2 * (o + p)); };

  for (size_t o = 0; o < ConvergenceOrder; o++) {
    for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
      for (size_t q = 0; q < misc::NumQuantities; q++) {
        qIPlus[o][q][p] = qP(o, q, p);
        qIMinus[o][q][p] = qM(o, q, p);
      }
      tractionResults.traction1[o][p] = t1(o, p);
      tractionResults.traction2[o][p] = t2(o, p);
    }
  }

  SUBCASE("Precompute Stress") {
    friction_law::common::precomputeStressFromQInterpolated(
        faultStresses, impAndEta, impMats, qInterpolatedPlus, qInterpolatedMinus);

    // Assure that qInterpolatedPlus and qInterpolatedMinus are const.
    for (size_t o = 0; o < ConvergenceOrder; o++) {
      for (size_t q = 0; q < misc::NumQuantities; q++) {
        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          REQUIRE(qIPlus[o][q][p] == qP(o, q, p));
          REQUIRE(qIMinus[o][q][p] == qM(o, q, p));
        }
      }
    }

    // Assure that faultstresses were computed correctly
    for (size_t o = 0; o < ConvergenceOrder; o++) {
      for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
        const real expectedNormalStress =
            impAndEta.etaP * (qM(o, 6, p) - qP(o, 6, p) + impAndEta.invZp * qP(o, 0, p) +
                              impAndEta.invZpNeig * qM(o, 0, p));
        const real expectedTraction1 =
            impAndEta.etaS * (qM(o, 7, p) - qP(o, 7, p) + impAndEta.invZs * qP(o, 3, p) +
                              impAndEta.invZsNeig * qM(o, 3, p));
        const real expectedTraction2 =
            impAndEta.etaS * (qM(o, 8, p) - qP(o, 8, p) + impAndEta.invZs * qP(o, 5, p) +
                              impAndEta.invZsNeig * qM(o, 5, p));
        REQUIRE(faultStresses.normalStress[o][p] ==
                AbsApprox(expectedNormalStress).epsilon(Epsilon));
        REQUIRE(faultStresses.traction1[o][p] == AbsApprox(expectedTraction1).epsilon(Epsilon));
        REQUIRE(faultStresses.traction2[o][p] == AbsApprox(expectedTraction2).epsilon(Epsilon));
      }
    }
  }

  SUBCASE("Postcompute Imposed State") {
    friction_law::common::postcomputeImposedStateFromNewStress(faultStresses,
                                                               tractionResults,
                                                               impAndEta,
                                                               impMats,
                                                               imposedStatePlus,
                                                               imposedStateMinus,
                                                               qInterpolatedPlus,
                                                               qInterpolatedMinus,
                                                               timeWeights);
    for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
      // index 0: Minus side
      // index 1: Plus side
      real expectedNormalStress[2]{};
      real expectedTraction1[2]{};
      real expectedTraction2[2]{};
      real expectedU[2]{};
      real expectedV[2]{};
      real expectedW[2]{};
      for (size_t o = 0; o < ConvergenceOrder; o++) {
        expectedNormalStress[0] += timeWeights[o] * faultStresses.normalStress[o][p];
        expectedTraction1[0] += timeWeights[o] * t1(o, p);
        expectedTraction2[0] += timeWeights[o] * t2(o, p);
        expectedU[0] +=
            timeWeights[o] *
            (qM(o, 6, p) - impAndEta.invZpNeig * (faultStresses.normalStress[o][p] - qM(o, 0, p)));
        expectedV[0] +=
            timeWeights[o] * (qM(o, 7, p) - impAndEta.invZsNeig * (t1(o, p) - qM(o, 3, p)));
        expectedW[0] +=
            timeWeights[o] * (qM(o, 8, p) - impAndEta.invZsNeig * (t2(o, p) - qM(o, 5, p)));

        expectedNormalStress[1] += timeWeights[o] * faultStresses.normalStress[o][p];
        expectedTraction1[1] += timeWeights[o] * t1(o, p);
        expectedTraction2[1] += timeWeights[o] * t2(o, p);
        expectedU[1] +=
            timeWeights[o] *
            (qP(o, 6, p) + impAndEta.invZp * (faultStresses.normalStress[o][p] - qP(o, 0, p)));
        expectedV[1] += timeWeights[o] * (qP(o, 7, p) + impAndEta.invZs * (t1(o, p) - qP(o, 3, p)));
        expectedW[1] += timeWeights[o] * (qP(o, 8, p) + impAndEta.invZs * (t2(o, p) - qP(o, 5, p)));
      }
      REQUIRE(iSMinus[0][p] == AbsApprox(expectedNormalStress[0]).epsilon(Epsilon));
      REQUIRE(iSMinus[3][p] == AbsApprox(expectedTraction1[0]).epsilon(Epsilon));
      REQUIRE(iSMinus[5][p] == AbsApprox(expectedTraction2[0]).epsilon(Epsilon));
      REQUIRE(iSMinus[6][p] == AbsApprox(expectedU[0]).epsilon(Epsilon));
      REQUIRE(iSMinus[7][p] == AbsApprox(expectedV[0]).epsilon(Epsilon));
      REQUIRE(iSMinus[8][p] == AbsApprox(expectedW[0]).epsilon(Epsilon));
      REQUIRE(iSPlus[0][p] == AbsApprox(expectedNormalStress[1]).epsilon(Epsilon));
      REQUIRE(iSPlus[3][p] == AbsApprox(expectedTraction1[1]).epsilon(Epsilon));
      REQUIRE(iSPlus[5][p] == AbsApprox(expectedTraction2[1]).epsilon(Epsilon));
      REQUIRE(iSPlus[6][p] == AbsApprox(expectedU[1]).epsilon(Epsilon));
      REQUIRE(iSPlus[7][p] == AbsApprox(expectedV[1]).epsilon(Epsilon));
      REQUIRE(iSPlus[8][p] == AbsApprox(expectedW[1]).epsilon(Epsilon));
    }
  }
}

} // namespace seissol::unit_test

#endif

#ifndef SEISSOL_FRICTIONSOLVERCOMMON_T_H
#define SEISSOL_FRICTIONSOLVERCOMMON_T_H

#include <numeric>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"

namespace seissol::unit_test::dr {

using namespace seissol;
using namespace seissol::dr;

TEST_CASE("Friction Solver Common") {
  FaultStresses faultStresses{};
  TractionResults tractionResults{};
  ImpedancesAndEta impAndEta;
  alignas(ALIGNMENT) real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  alignas(ALIGNMENT) real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  // TODO(NONLINEAR) Reset this file to the state of ff11d9426539be5a3cbb23787756b7f4050d5cc3

  alignas(ALIGNMENT) real qStrainInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  alignas(ALIGNMENT)
      real qStrainInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  alignas(ALIGNMENT) real imposedStatePlus[tensor::QInterpolated::size()];
  alignas(ALIGNMENT) real imposedStateMinus[tensor::QInterpolated::size()];
  double timeWeights[CONVERGENCE_ORDER];
  std::iota(std::begin(timeWeights), std::end(timeWeights), 1);

  using QInterpolatedShapeT = real(*)[misc::numQuantities][misc::numPaddedPoints];
  auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
  auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));
  using ImposedStateShapeT = real(*)[misc::numPaddedPoints];
  auto* iSPlus = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
  auto* iSMinus = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

  impAndEta.etaP = 2.0;
  impAndEta.etaS = 4.0;
  impAndEta.zp = 10.0;
  impAndEta.zs = 20.0;
  impAndEta.zpNeig = 15.0;
  impAndEta.zsNeig = 25.0;
  impAndEta.invZp = 1.0 / impAndEta.zp;
  impAndEta.invZs = 1.0 / impAndEta.zs;
  impAndEta.invZpNeig = 1.0 / impAndEta.zpNeig;
  impAndEta.invZsNeig = 1.0 / impAndEta.zsNeig;

  auto qP = [](size_t o, size_t q, size_t p) { return static_cast<real>(o + q + p); };
  auto qM = [](size_t o, size_t q, size_t p) { return static_cast<real>(2 * (o + q + p)); };
  auto t1 = [](size_t o, size_t p) { return static_cast<real>(o + p); };
  auto t2 = [](size_t o, size_t p) { return static_cast<real>(2 * (o + p)); };

  for (size_t o = 0; o < CONVERGENCE_ORDER; o++) {
    for (size_t p = 0; p < misc::numPaddedPoints; p++) {
      for (size_t q = 0; q < 9; q++) {
        qIPlus[o][q][p] = qP(o, q, p);
        qIMinus[o][q][p] = qM(o, q, p);
      }
      tractionResults.traction1[o][p] = t1(o, p);
      tractionResults.traction2[o][p] = t2(o, p);
    }
  }

  SUBCASE("Precompute Stress") {
    friction_law::common::precomputeStressFromQInterpolated<
        seissol::dr::friction_law::common::RangeType::CPU>(faultStresses,
                                                           impAndEta,
                                                           qInterpolatedPlus,
                                                           qInterpolatedMinus,
                                                           qStrainInterpolatedMinus,
                                                           qStrainInterpolatedPlus);

    // Assure that qInterpolatedPlus and qInterpolatedMinus are const.
    for (size_t o = 0; o < CONVERGENCE_ORDER; o++) {
      for (size_t q = 0; q < 9; q++) {
        for (size_t p = 0; p < misc::numPaddedPoints; p++) {
          REQUIRE(qIPlus[o][q][p] == qP(o, q, p));
          REQUIRE(qIMinus[o][q][p] == qM(o, q, p));
        }
      }
    }

    // Assure that faultstresses were computed correctly
    for (size_t o = 0; o < CONVERGENCE_ORDER; o++) {
      for (size_t p = 0; p < misc::numPaddedPoints; p++) {
        real expectedNormalStress =
            impAndEta.etaP * (qM(o, 6, p) - qP(o, 6, p) + impAndEta.invZp * qP(o, 0, p) +
                              impAndEta.invZpNeig * qM(o, 0, p));
        real expectedTraction1 =
            impAndEta.etaS * (qM(o, 7, p) - qP(o, 7, p) + impAndEta.invZs * qP(o, 3, p) +
                              impAndEta.invZsNeig * qM(o, 3, p));
        real expectedTraction2 =
            impAndEta.etaS * (qM(o, 8, p) - qP(o, 8, p) + impAndEta.invZs * qP(o, 5, p) +
                              impAndEta.invZsNeig * qM(o, 5, p));
        REQUIRE(faultStresses.normalStress[o][p] == expectedNormalStress);
        REQUIRE(faultStresses.traction1[o][p] == expectedTraction1);
        REQUIRE(faultStresses.traction2[o][p] == expectedTraction2);
      }
    }
  }

  SUBCASE("Postcompute Imposed State") {
    friction_law::common::postcomputeImposedStateFromNewStress(faultStresses,
                                                               tractionResults,
                                                               impAndEta,
                                                               imposedStatePlus,
                                                               imposedStateMinus,
                                                               qInterpolatedPlus,
                                                               qInterpolatedMinus,
                                                               timeWeights);
    for (size_t p = 0; p < misc::numPaddedPoints; p++) {
      // index 0: Minus side
      // index 1: Plus side
      real expectedNormalStress[2]{};
      real expectedTraction1[2]{};
      real expectedTraction2[2]{};
      real expectedU[2]{};
      real expectedV[2]{};
      real expectedW[2]{};
      for (size_t o = 0; o < CONVERGENCE_ORDER; o++) {
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
      REQUIRE(iSMinus[0][p] == expectedNormalStress[0]);
      REQUIRE(iSMinus[3][p] == expectedTraction1[0]);
      REQUIRE(iSMinus[5][p] == expectedTraction2[0]);
      REQUIRE(iSMinus[6][p] == expectedU[0]);
      REQUIRE(iSMinus[7][p] == expectedV[0]);
      REQUIRE(iSMinus[8][p] == expectedW[0]);
      REQUIRE(iSPlus[0][p] == expectedNormalStress[1]);
      REQUIRE(iSPlus[3][p] == expectedTraction1[1]);
      REQUIRE(iSPlus[5][p] == expectedTraction2[1]);
      REQUIRE(iSPlus[6][p] == expectedU[1]);
      REQUIRE(iSPlus[7][p] == expectedV[1]);
      REQUIRE(iSPlus[8][p] == expectedW[1]);
    }
  }
}

} // namespace seissol::unit_test::dr

#endif // SEISSOL_FRICTIONSOLVERCOMMON_T_H

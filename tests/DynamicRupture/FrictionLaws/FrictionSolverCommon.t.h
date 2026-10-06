// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "TestHelper.h"

#include <cstdint>
#include <numeric>

namespace seissol::unit_test {

using namespace seissol;
using namespace seissol::dr;

TEST_CASE("Friction Solver Common" * doctest::test_suite("dynamicrupture")) {
  if constexpr (model::MaterialT::SupportsDR) {
    FaultStresses<Executor::Host> faultStresses{};
    TractionResults<Executor::Host> tractionResults{};
    ImposedState<Executor::Host> imposedState{};
    ImpedancesAndEta impAndEta;
    alignas(Alignment)
        real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()] = {{}};
    alignas(Alignment)
        real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()] = {{}};
    alignas(Alignment) real imposedStatePlus[tensor::QInterpolated::size()] = {};
    alignas(Alignment) real imposedStateMinus[tensor::QInterpolated::size()] = {};
    real timeWeights[misc::TimeSteps]{};
    std::iota(std::begin(timeWeights), std::end(timeWeights), 1);
    constexpr real Epsilon = 1e6 * std::numeric_limits<real>::epsilon();

    constexpr auto GpuRange = friction_law::common::RangeType::GPU;

    using QInterpolatedShapeT = real(*)[misc::NumQuantities][misc::NumPaddedPoints];
    auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus));
    auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus));
    using ImposedStateShapeT = real(*)[misc::NumPaddedPoints];
    auto* iSPlus = reinterpret_cast<ImposedStateShapeT>(imposedStatePlus);
    auto* iSMinus = reinterpret_cast<ImposedStateShapeT>(imposedStateMinus);

    const real zp = 10.0;
    const real zs = 20.0;
    const real zpNeig = 15.0;
    const real zsNeig = 25.0;
    const real etaP = zp * zpNeig / (zp + zpNeig);
    const real etaS = zs * zsNeig / (zs + zsNeig);
    const real invZp = 1.0 / zp;
    const real invZs = 1.0 / zs;
    const real invZpNeig = 1.0 / zpNeig;
    const real invZsNeig = 1.0 / zsNeig;
    impAndEta.zp.fill(zp);
    impAndEta.zs.fill(zs);
    impAndEta.zpNeig.fill(zpNeig);
    impAndEta.zsNeig.fill(zsNeig);
    impAndEta.etaP.fill(etaP);
    impAndEta.etaS.fill(etaS);
    impAndEta.invZp.fill(invZp);
    impAndEta.invZs.fill(invZs);
    impAndEta.invZpNeig.fill(invZpNeig);
    impAndEta.invZsNeig.fill(invZsNeig);

    ImpedanceMatrices impMats;
    for (std::size_t point = 0; point < ImpedancePoints; ++point) {
      auto etaView = init::eta::view::create(impMats.eta.at(point));
      etaView(0, 0) = etaP;
      etaView(1, 1) = etaS;
      etaView(2, 2) = etaS;
      auto impedanceView = init::Zplus::view::create(impMats.impedance.at(point));
      impedanceView(0, 0) = invZp;
      impedanceView(1, 1) = invZs;
      impedanceView(2, 2) = invZs;
      auto impedanceNeigView = init::Zminus::view::create(impMats.impedanceNeig.at(point));
      impedanceNeigView(0, 0) = invZpNeig;
      impedanceNeigView(1, 1) = invZsNeig;
      impedanceNeigView(2, 2) = invZsNeig;
    }

    auto qP = [](size_t o, size_t q, size_t p) { return static_cast<real>(o + q + p); };
    auto qM = [](size_t o, size_t q, size_t p) { return static_cast<real>(2 * (o + q + p)); };
    auto t1 = [](size_t o, size_t p) { return static_cast<real>(o + p); };
    auto t2 = [](size_t o, size_t p) { return static_cast<real>(2 * (o + p)); };
    // deliberately different from faultStresses.normalStress, so that the postcompute test
    // discriminates between the trial and the post-solve normal stress
    auto tn = [](size_t o, size_t p) { return static_cast<real>(3 * (o + p) + 1); };

    // qInterpolated still carries every time step; FaultStresses and TractionResults do not --
    // they hold a single time slice and are refilled on every step of the friction loop.
    for (size_t o = 0; o < misc::TimeSteps; o++) {
      for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
        for (size_t q = 0; q < misc::NumQuantities; q++) {
          qIPlus[o][q][p] = qP(o, q, p);
          qIMinus[o][q][p] = qM(o, q, p);
        }
      }
    }

    // stands in for the friction law: fills the single-slice traction results for step o
    auto seedTractionResults = [&](size_t o) {
      for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
        tractionResults.normalStress[p] = tn(o, p);
        tractionResults.traction1[p] = t1(o, p);
        tractionResults.traction2[p] = t2(o, p);
      }
    };

    // the trial stresses precomputeStressFromQInterpolated is expected to produce for step o
    auto trialNormalStress = [&](size_t o, size_t p) {
      return etaP * (qM(o, 6, p) - qP(o, 6, p) + invZp * qP(o, 0, p) + invZpNeig * qM(o, 0, p));
    };
    auto trialTraction1 = [&](size_t o, size_t p) {
      return etaS * (qM(o, 7, p) - qP(o, 7, p) + invZs * qP(o, 3, p) + invZsNeig * qM(o, 3, p));
    };
    auto trialTraction2 = [&](size_t o, size_t p) {
      return etaS * (qM(o, 8, p) - qP(o, 8, p) + invZs * qP(o, 5, p) + invZsNeig * qM(o, 5, p));
    };

    SUBCASE("Precompute Stress") {
      for (size_t o = 0; o < misc::TimeSteps; o++) {
        friction_law::common::precomputeStressFromQInterpolated(
            faultStresses, impAndEta, impMats, qInterpolatedPlus, qInterpolatedMinus, 1.0, o);

        // Assure that the faultstresses of *this* step were computed correctly. Since the struct
        // holds a single slice, a step index that is ignored somewhere would show up right here.
        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          CHECK(faultStresses.normalStress[p] ==
                AbsApprox(trialNormalStress(o, p)).epsilon(Epsilon));
          CHECK(faultStresses.traction1[p] == AbsApprox(trialTraction1(o, p)).epsilon(Epsilon));
          CHECK(faultStresses.traction2[p] == AbsApprox(trialTraction2(o, p)).epsilon(Epsilon));
        }
      }

      // Assure that qInterpolatedPlus and qInterpolatedMinus are const.
      for (size_t o = 0; o < misc::TimeSteps; o++) {
        for (size_t q = 0; q < misc::NumQuantities; q++) {
          for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
            CHECK(qIPlus[o][q][p] == qP(o, q, p));
            CHECK(qIMinus[o][q][p] == qM(o, q, p));
          }
        }
      }
    }

    SUBCASE("Initialize Traction Results") {
      // the seeding has to happen on every step now, not once up front -- comparing against the
      // trial stress of the current step catches a seed that is left over from an earlier one
      for (size_t o = 0; o < misc::TimeSteps; o++) {
        seedTractionResults(o);
        friction_law::common::precomputeStressFromQInterpolated(
            faultStresses, impAndEta, impMats, qInterpolatedPlus, qInterpolatedMinus, 1.0, o);
        friction_law::common::initializeTractionResults(faultStresses, tractionResults);

        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          // the trial normal stress of *this* step is seeded ...
          CHECK(tractionResults.normalStress[p] ==
                AbsApprox(trialNormalStress(o, p)).epsilon(Epsilon));
          // ... i.e. the value put there beforehand is really gone ...
          CHECK(tractionResults.normalStress[p] != AbsApprox(tn(o, p)).epsilon(Epsilon));
          // ... and nothing else is touched
          CHECK(tractionResults.traction1[p] == AbsApprox(t1(o, p)).epsilon(Epsilon));
          CHECK(tractionResults.traction2[p] == AbsApprox(t2(o, p)).epsilon(Epsilon));
        }
      }
    }

    SUBCASE("Postcompute Imposed State") {
      // finalizeImposedState has to overwrite the output, not accumulate into it
      for (size_t i = 0; i < tensor::QInterpolated::size(); i++) {
        imposedStatePlus[i] = static_cast<real>(-1.0);
        imposedStateMinus[i] = static_cast<real>(-1.0);
      }

      for (size_t o = 0; o < misc::TimeSteps; o++) {
        seedTractionResults(o);
        friction_law::common::postcomputeImposedStateFromNewStress(imposedState,
                                                                   faultStresses,
                                                                   tractionResults,
                                                                   impAndEta,
                                                                   impMats,
                                                                   qInterpolatedPlus,
                                                                   qInterpolatedMinus,
                                                                   o,
                                                                   timeWeights[o]);
      }
      friction_law::common::finalizeImposedState(imposedState, imposedStatePlus, imposedStateMinus);

      for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
        // index 0: Minus side
        // index 1: Plus side
        real expectedNormalStress[2]{};
        real expectedTraction1[2]{};
        real expectedTraction2[2]{};
        real expectedU[2]{};
        real expectedV[2]{};
        real expectedW[2]{};
        for (size_t o = 0; o < misc::TimeSteps; o++) {
          expectedNormalStress[0] += timeWeights[o] * tn(o, p);
          expectedTraction1[0] += timeWeights[o] * t1(o, p);
          expectedTraction2[0] += timeWeights[o] * t2(o, p);
          expectedU[0] += timeWeights[o] * (qM(o, 6, p) - invZpNeig * (tn(o, p) - qM(o, 0, p)));
          expectedV[0] += timeWeights[o] * (qM(o, 7, p) - invZsNeig * (t1(o, p) - qM(o, 3, p)));
          expectedW[0] += timeWeights[o] * (qM(o, 8, p) - invZsNeig * (t2(o, p) - qM(o, 5, p)));

          expectedNormalStress[1] += timeWeights[o] * tn(o, p);
          expectedTraction1[1] += timeWeights[o] * t1(o, p);
          expectedTraction2[1] += timeWeights[o] * t2(o, p);
          expectedU[1] += timeWeights[o] * (qP(o, 6, p) + invZp * (tn(o, p) - qP(o, 0, p)));
          expectedV[1] += timeWeights[o] * (qP(o, 7, p) + invZs * (t1(o, p) - qP(o, 3, p)));
          expectedW[1] += timeWeights[o] * (qP(o, 8, p) + invZs * (t2(o, p) - qP(o, 5, p)));
        }
        CHECK(iSMinus[0][p] == AbsApprox(expectedNormalStress[0]).epsilon(Epsilon));
        CHECK(iSMinus[3][p] == AbsApprox(expectedTraction1[0]).epsilon(Epsilon));
        CHECK(iSMinus[5][p] == AbsApprox(expectedTraction2[0]).epsilon(Epsilon));
        CHECK(iSMinus[6][p] == AbsApprox(expectedU[0]).epsilon(Epsilon));
        CHECK(iSMinus[7][p] == AbsApprox(expectedV[0]).epsilon(Epsilon));
        CHECK(iSMinus[8][p] == AbsApprox(expectedW[0]).epsilon(Epsilon));
        CHECK(iSPlus[0][p] == AbsApprox(expectedNormalStress[1]).epsilon(Epsilon));
        CHECK(iSPlus[3][p] == AbsApprox(expectedTraction1[1]).epsilon(Epsilon));
        CHECK(iSPlus[5][p] == AbsApprox(expectedTraction2[1]).epsilon(Epsilon));
        CHECK(iSPlus[6][p] == AbsApprox(expectedU[1]).epsilon(Epsilon));
        CHECK(iSPlus[7][p] == AbsApprox(expectedV[1]).epsilon(Epsilon));
        CHECK(iSPlus[8][p] == AbsApprox(expectedW[1]).epsilon(Epsilon));

        // YY, ZZ, YZ take no part in the fault-normal Riemann problem; the accumulator never
        // touches them, so finalizeImposedState has to write the zeros through
        CHECK(iSMinus[1][p] == AbsApprox(0.0).epsilon(Epsilon));
        CHECK(iSMinus[2][p] == AbsApprox(0.0).epsilon(Epsilon));
        CHECK(iSMinus[4][p] == AbsApprox(0.0).epsilon(Epsilon));
        CHECK(iSPlus[1][p] == AbsApprox(0.0).epsilon(Epsilon));
        CHECK(iSPlus[2][p] == AbsApprox(0.0).epsilon(Epsilon));
        CHECK(iSPlus[4][p] == AbsApprox(0.0).epsilon(Epsilon));
      }
    }

    SUBCASE("Subtract Own Traces") {
      // the jump against the own trace of either side, the trace integrated with the time weights
      ImposedState<Executor::Host> state{};
      auto before = [](int side, size_t q, size_t p) {
        return static_cast<real>(side * (1.0 + q) + 0.01 * p);
      };
      for (size_t q = 0; q < misc::NumQuantities; q++) {
        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          state.plus[q][p] = before(1, q, p);
          state.minus[q][p] = before(-2, q, p);
        }
      }
      friction_law::common::subtractOwnTraces(
          state, qInterpolatedPlus, qInterpolatedMinus, timeWeights);
      for (size_t q = 0; q < misc::NumQuantities; q++) {
        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          real tracePlus = 0;
          real traceMinus = 0;
          for (size_t o = 0; o < misc::TimeSteps; o++) {
            tracePlus += timeWeights[o] * qP(o, q, p);
            traceMinus += timeWeights[o] * qM(o, q, p);
          }
          CHECK(state.plus[q][p] == AbsApprox(before(1, q, p) - tracePlus).epsilon(Epsilon));
          CHECK(state.minus[q][p] == AbsApprox(before(-2, q, p) - traceMinus).epsilon(Epsilon));
        }
      }
    }

    SUBCASE("Device Range Matches Host Range") {
      // The device specializations collapse every point-indexed array to a scalar and address the
      // padded point through startIndex instead. Instantiating them for RangeType::GPU on the host
      // is the only coverage that path gets in a CPU build, and it pins down the host/device index
      // handling that the single-slice rework touches.
      alignas(Alignment) real deviceImposedStatePlus[tensor::QInterpolated::size()] = {};
      alignas(Alignment) real deviceImposedStateMinus[tensor::QInterpolated::size()] = {};
      auto* dSPlus = reinterpret_cast<ImposedStateShapeT>(deviceImposedStatePlus);
      auto* dSMinus = reinterpret_cast<ImposedStateShapeT>(deviceImposedStateMinus);

      for (size_t o = 0; o < misc::TimeSteps; o++) {
        friction_law::common::precomputeStressFromQInterpolated(
            faultStresses, impAndEta, impMats, qInterpolatedPlus, qInterpolatedMinus, 1.0, o);
        friction_law::common::initializeTractionResults(faultStresses, tractionResults);
        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          tractionResults.traction1[p] = t1(o, p);
          tractionResults.traction2[p] = t2(o, p);
        }
        friction_law::common::postcomputeImposedStateFromNewStress(imposedState,
                                                                   faultStresses,
                                                                   tractionResults,
                                                                   impAndEta,
                                                                   impMats,
                                                                   qInterpolatedPlus,
                                                                   qInterpolatedMinus,
                                                                   o,
                                                                   timeWeights[o]);
      }
      friction_law::common::finalizeImposedState(imposedState, imposedStatePlus, imposedStateMinus);

      for (std::uint32_t p = 0; p < misc::NumPaddedPoints; p++) {
        FaultStresses<Executor::Device> deviceFaultStresses{};
        TractionResults<Executor::Device> deviceTractionResults{};
        ImposedState<Executor::Device> deviceImposedState{};

        for (std::uint32_t o = 0; o < misc::TimeSteps; o++) {
          friction_law::common::precomputeStressFromQInterpolated<GpuRange>(deviceFaultStresses,
                                                                            impAndEta,
                                                                            impMats,
                                                                            qInterpolatedPlus,
                                                                            qInterpolatedMinus,
                                                                            1.0,
                                                                            o,
                                                                            p);
          friction_law::common::initializeTractionResults<GpuRange>(
              deviceFaultStresses, deviceTractionResults, p);
          deviceTractionResults.traction1 = t1(o, p);
          deviceTractionResults.traction2 = t2(o, p);
          friction_law::common::postcomputeImposedStateFromNewStress<GpuRange>(
              deviceImposedState,
              deviceFaultStresses,
              deviceTractionResults,
              impAndEta,
              impMats,
              qInterpolatedPlus,
              qInterpolatedMinus,
              o,
              timeWeights[o],
              p);
        }
        friction_law::common::finalizeImposedState<GpuRange>(
            deviceImposedState, deviceImposedStatePlus, deviceImposedStateMinus, p);
      }

      for (size_t q = 0; q < misc::NumQuantities; q++) {
        for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
          CHECK(dSPlus[q][p] == AbsApprox(iSPlus[q][p]).epsilon(Epsilon).delta(Epsilon));
          CHECK(dSMinus[q][p] == AbsApprox(iSMinus[q][p]).epsilon(Epsilon).delta(Epsilon));
        }
      }
    }
  }
}

} // namespace seissol::unit_test

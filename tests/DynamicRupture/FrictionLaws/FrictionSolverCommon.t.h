// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "Numerical/Quadrature.h"
#include "TestHelper.h"

#include <algorithm>
#include <cstdint>
#include <limits>
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
      return impAndEta.etaP * (qM(o, 6, p) - qP(o, 6, p) + impAndEta.invZp * qP(o, 0, p) +
                               impAndEta.invZpNeig * qM(o, 0, p));
    };
    auto trialTraction1 = [&](size_t o, size_t p) {
      return impAndEta.etaS * (qM(o, 7, p) - qP(o, 7, p) + impAndEta.invZs * qP(o, 3, p) +
                               impAndEta.invZsNeig * qM(o, 3, p));
    };
    auto trialTraction2 = [&](size_t o, size_t p) {
      return impAndEta.etaS * (qM(o, 8, p) - qP(o, 8, p) + impAndEta.invZs * qP(o, 5, p) +
                               impAndEta.invZsNeig * qM(o, 5, p));
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
          expectedU[0] +=
              timeWeights[o] * (qM(o, 6, p) - impAndEta.invZpNeig * (tn(o, p) - qM(o, 0, p)));
          expectedV[0] +=
              timeWeights[o] * (qM(o, 7, p) - impAndEta.invZsNeig * (t1(o, p) - qM(o, 3, p)));
          expectedW[0] +=
              timeWeights[o] * (qM(o, 8, p) - impAndEta.invZsNeig * (t2(o, p) - qM(o, 5, p)));

          expectedNormalStress[1] += timeWeights[o] * tn(o, p);
          expectedTraction1[1] += timeWeights[o] * t1(o, p);
          expectedTraction2[1] += timeWeights[o] * t2(o, p);
          expectedU[1] +=
              timeWeights[o] * (qP(o, 6, p) + impAndEta.invZp * (tn(o, p) - qP(o, 0, p)));
          expectedV[1] +=
              timeWeights[o] * (qP(o, 7, p) + impAndEta.invZs * (t1(o, p) - qP(o, 3, p)));
          expectedW[1] +=
              timeWeights[o] * (qP(o, 8, p) + impAndEta.invZs * (t2(o, p) - qP(o, 5, p)));
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

    SUBCASE("Friction energy accumulates in double") {
      // late in a simulation, the sums of the energy output are large and their increments tiny;
      // in single precision, the increments would be rounded away
      DREnergyOutput energyOutput{};
      real slipRateMagnitude[misc::NumPaddedPoints]{};
      for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
        energyOutput.accumulatedSlip[p] = 1.0;
        slipRateMagnitude[p] = static_cast<real>(1e-9);
      }
      real spaceWeights[seissol::kernels::NumSpaceQuadraturePoints]{};
      const DRGodunovData godunovData{};
      friction_law::common::computeFrictionEnergy(energyOutput,
                                                  qInterpolatedPlus,
                                                  qInterpolatedMinus,
                                                  impAndEta,
                                                  timeWeights,
                                                  spaceWeights,
                                                  godunovData,
                                                  slipRateMagnitude,
                                                  false);
      double gain = 0;
      for (size_t o = 0; o < misc::TimeSteps; o++) {
        gain += static_cast<double>(timeWeights[o] * slipRateMagnitude[0]);
      }
      for (size_t p = 0; p < misc::NumPaddedPoints; p++) {
        CHECK(energyOutput.accumulatedSlip[p] - 1.0 == doctest::Approx(gain).epsilon(1e-6));
      }
    }
  }
}

TEST_CASE("Friction solver event times" * doctest::test_suite("dynamicrupture")) {
  using friction_law::common::crossingTime;
  const auto real0 = [](double x) { return static_cast<real>(x); };

  SUBCASE("A crossing is interpolated linearly, rising or falling") {
    CHECK(crossingTime(real0(0.0), real0(2.0), real0(1.0), real0(10.0), real0(12.0)) ==
          doctest::Approx(11.0));
    CHECK(crossingTime(real0(0.8), real0(0.4), real0(0.6), real0(1.0), real0(2.0)) ==
          doctest::Approx(1.5));
    // landing on the threshold exactly, the crossing is at the current time
    CHECK(crossingTime(real0(0.0), real0(1.0), real0(1.0), real0(1.0), real0(2.0)) ==
          doctest::Approx(2.0));
  }

  SUBCASE("A previous value past the threshold gives the previous time") {
    // moving away from the threshold, back toward it, or not at all
    CHECK(crossingTime(real0(2.0), real0(3.0), real0(1.0), real0(5.0), real0(6.0)) == 5.0);
    CHECK(crossingTime(real0(3.0), real0(2.0), real0(1.0), real0(5.0), real0(6.0)) == 5.0);
    CHECK(crossingTime(real0(2.0), real0(2.0), real0(1.0), real0(5.0), real0(6.0)) == 5.0);
    // and so continuously as the current value passes the previous one
    CHECK(crossingTime(real0(2.0), real0(2.001), real0(1.0), real0(5.0), real0(6.0)) == 5.0);
    CHECK(crossingTime(real0(2.0), real0(1.999), real0(1.0), real0(5.0), real0(6.0)) == 5.0);
  }

  // the time points of a step as the time stepping gets them
  const double stepLength = 0.7;
  const auto timePoints =
      seissol::quadrature::ShiftedGaussLegendre(misc::TimeSteps, 0, stepLength).first;
  const auto frictionTime = friction_law::FrictionSolver::computeDeltaT(timePoints);
  real deltaT[misc::TimeSteps]{};
  for (std::size_t i = 0; i < misc::TimeSteps; ++i) {
    deltaT[i] = static_cast<real>(frictionTime.deltaT[i]);
  }
  const double tolerance = 64 * std::numeric_limits<real>::epsilon();

  SUBCASE("Sub-step times follow the time points of the step") {
    const double stepStart = 3.0;
    // the last time point of a previous step of the same length
    const double previousStepPoint = stepStart - (stepLength - timePoints.back());
    for (uint32_t k = 0; k < misc::TimeSteps; ++k) {
      CAPTURE(k);
      const auto times = friction_law::common::subStepTimes(real0(stepStart), deltaT, k);
      CHECK(times.point == doctest::Approx(stepStart + timePoints[k]).epsilon(tolerance));
      const double previous = k == 0 ? previousStepPoint : stepStart + timePoints[k - 1];
      CHECK(times.previousPoint == doctest::Approx(previous).epsilon(tolerance));
      const double start = k == 0 ? stepStart : stepStart + timePoints[k - 1];
      CHECK(times.start == doctest::Approx(start).epsilon(tolerance));
      const double end =
          k + 1 < misc::TimeSteps ? stepStart + timePoints[k] : stepStart + stepLength;
      CHECK(times.end == doctest::Approx(end).epsilon(tolerance));
    }
  }

  SUBCASE("The first evaluation of a simulation has no previous one") {
    const auto first = friction_law::common::subStepTimes(real0(0.0), deltaT, 0);
    CHECK(first.previousPoint == 0.0);
    CHECK(first.point == 0.0);
    const auto second = friction_law::common::subStepTimes(real0(0.0), deltaT, 1);
    CHECK(second.previousPoint == doctest::Approx(timePoints[0]).epsilon(tolerance));

    // a prescribed slip rate keeps them too
    const auto rate = friction_law::common::slipRateTimes(first, true);
    CHECK(rate.previousPoint == 0.0);
    CHECK(rate.point == 0.0);
  }

  SUBCASE("A prescribed slip rate belongs to the ends of the sub-step interval") {
    for (uint32_t k = 0; k < misc::TimeSteps; ++k) {
      CAPTURE(k);
      const auto times = friction_law::common::subStepTimes(real0(3.0), deltaT, k);
      const auto rate = friction_law::common::slipRateTimes(times, true);
      CHECK(rate.previousPoint == times.start);
      CHECK(rate.point == times.end);
      const auto solved = friction_law::common::slipRateTimes(times, false);
      CHECK(solved.previousPoint == times.previousPoint);
      CHECK(solved.point == times.point);
    }
  }

  // per point: a previous slip rate below the threshold of 1 mm/s, and a current one above it
  // for all points but the last
  bool pending[misc::NumPaddedPoints]{};
  real ruptureTime[misc::NumPaddedPoints]{};
  real previousSlipRate[misc::NumPaddedPoints]{};
  real slipRate[misc::NumPaddedPoints]{};
  for (std::size_t p = 0; p < misc::NumPaddedPoints; ++p) {
    pending[p] = true;
    previousSlipRate[p] = static_cast<real>(0.0005 * p / misc::NumPaddedPoints);
    slipRate[p] = p + 1 < misc::NumPaddedPoints ? static_cast<real>(0.001 * (2 + p)) : 0.0;
  }

  SUBCASE("The rupture time moves with the crossing") {
    friction_law::common::saveRuptureFrontOutput(
        pending, ruptureTime, previousSlipRate, slipRate, real0(2.0), real0(3.0));
    for (std::size_t p = 0; p + 1 < misc::NumPaddedPoints; ++p) {
      CAPTURE(p);
      CHECK_FALSE(pending[p]);
      CHECK(ruptureTime[p] ==
            crossingTime(previousSlipRate[p], slipRate[p], real0(0.001), real0(2.0), real0(3.0)));
      CHECK(ruptureTime[p] > 2.0);
      CHECK(ruptureTime[p] <= 3.0);
    }
    CHECK(pending[misc::NumPaddedPoints - 1]);
    CHECK(ruptureTime[misc::NumPaddedPoints - 1] == 0.0);
  }

  SUBCASE("A rupture at the first evaluation of a simulation is at 0") {
    friction_law::common::saveRuptureFrontOutput(
        pending, ruptureTime, previousSlipRate, slipRate, real0(0.0), real0(0.0));
    CHECK_FALSE(pending[0]);
    CHECK(ruptureTime[0] == 0.0);
  }

  SUBCASE("The device range matches the host range") {
    // the device variant takes the previous slip rate of its point alone
    bool devicePending[misc::NumPaddedPoints]{};
    real deviceRuptureTime[misc::NumPaddedPoints]{};
    std::copy_n(pending, misc::NumPaddedPoints, devicePending);
    for (uint32_t p = 0; p < misc::NumPaddedPoints; ++p) {
      const real previous = previousSlipRate[p];
      friction_law::common::saveRuptureFrontOutput<friction_law::common::RangeType::GPU>(
          devicePending, deviceRuptureTime, &previous, slipRate, real0(2.0), real0(3.0), p);
    }
    friction_law::common::saveRuptureFrontOutput(
        pending, ruptureTime, previousSlipRate, slipRate, real0(2.0), real0(3.0));
    for (std::size_t p = 0; p < misc::NumPaddedPoints; ++p) {
      CAPTURE(p);
      CHECK(devicePending[p] == pending[p]);
      CHECK(deviceRuptureTime[p] == ruptureTime[p]);
    }
  }
}

} // namespace seissol::unit_test

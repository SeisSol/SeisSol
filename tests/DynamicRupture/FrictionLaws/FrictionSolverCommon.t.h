// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/FrictionLaws/FrictionSolverCommon.h"
#include "DynamicRupture/Misc.h"
#include "TestHelper.h"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>

namespace seissol::unit_test {

using namespace seissol;
using namespace seissol::dr;

TEST_CASE("Friction Solver Common" * doctest::test_suite("dynamicrupture")) {
  // Every expectation below is written for an impedance that belongs to the face: the scalars of
  // ImpedancesAndEta, and the matrices assembled from them once for all of its nodes. Where the
  // impedance follows the state of a node instead, none of them is the value to compare against.
  if constexpr (model::MaterialT::SupportsDR && !NodalImpedance) {
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

    // Zeroed, and never read: the precomputation takes this up only where the impedance follows
    // the state of a node, which the guard above excludes.
    NodalImpedanceParameters nodalImpedanceParams{};

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
        friction_law::common::precomputeStressFromQInterpolated(faultStresses,
                                                                impAndEta,
                                                                impMats,
                                                                &nodalImpedanceParams,
                                                                qInterpolatedPlus,
                                                                qInterpolatedMinus,
                                                                1.0,
                                                                o);

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
        friction_law::common::precomputeStressFromQInterpolated(faultStresses,
                                                                impAndEta,
                                                                impMats,
                                                                &nodalImpedanceParams,
                                                                qInterpolatedPlus,
                                                                qInterpolatedMinus,
                                                                1.0,
                                                                o);
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
        friction_law::common::precomputeStressFromQInterpolated(faultStresses,
                                                                impAndEta,
                                                                impMats,
                                                                &nodalImpedanceParams,
                                                                qInterpolatedPlus,
                                                                qInterpolatedMinus,
                                                                1.0,
                                                                o);
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
                                                                            &nodalImpedanceParams,
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

namespace frictionenergy {

/// The matrices of node p at step o. Nothing in particular -- only different from node to node
/// and from step to step, and not symmetric, so that a read of the wrong node, of the wrong step
/// or of the transpose all change the result.
inline void nodeMatrices(std::size_t o, std::size_t p, real plus[9], real minus[9], real eta[9]) {
  for (std::size_t k = 0; k < 9; ++k) {
    plus[k] = static_cast<real>(1.0e-3 * (1.0 + 0.7 * k + 0.05 * p + 0.2 * o));
    minus[k] = static_cast<real>(1.0e-3 * (2.0 + 0.3 * k * k + 0.02 * p + 0.1 * o));
    eta[k] = static_cast<real>(1.0e2 * (1.0 + 0.4 * k + 0.01 * k * p + 0.3 * o));
  }
}

/// Puts the matrices of step o where the friction solve keeps them: with the stresses of every
/// node, where the impedance belongs to the node. A template, so that the members it fills are
/// looked up once the stresses are put in: they have them only where it belongs to the node.
template <typename StressesT>
void fillHost([[maybe_unused]] StressesT& stresses, [[maybe_unused]] std::size_t o) {
  if constexpr (NodalImpedance) {
    for (std::size_t p = 0; p < misc::NumPaddedPoints; ++p) {
      real plus[9];
      real minus[9];
      real eta[9];
      nodeMatrices(o, p, plus, minus, eta);
      for (std::size_t k = 0; k < 9; ++k) {
        stresses.admittance[k][p] = plus[k];
        stresses.admittanceNeig[k][p] = minus[k];
        stresses.eta[k][p] = eta[k];
      }
    }
  }
}

template <typename StressesT>
void fillDevice([[maybe_unused]] StressesT& stresses,
                [[maybe_unused]] std::size_t o,
                [[maybe_unused]] std::size_t p) {
  if constexpr (NodalImpedance) {
    real plus[9];
    real minus[9];
    real eta[9];
    nodeMatrices(o, p, plus, minus, eta);
    for (std::size_t k = 0; k < 9; ++k) {
      stresses.admittance[k] = plus[k];
      stresses.admittanceNeig[k] = minus[k];
      stresses.eta[k] = eta[k];
    }
  }
}

} // namespace frictionenergy

TEST_CASE("Frictional energy weighs the two traces the way the Riemann problem does" *
          doctest::skip(!model::MaterialT::SupportsDR ||
                        model::MaterialT::Type == model::MaterialType::Anisotropic) *
          doctest::test_suite("dynamicrupture")) {
  // An anisotropic face weighs with the matrices of the face, as assembled for its averaging
  // kernel; that path is not the subject here.
  using namespace frictionenergy;
  using namespace misc::quantity_indices;
  namespace common = friction_law::common;
  constexpr auto GpuRange = common::RangeType::GPU;

  alignas(Alignment) real qInterpolatedPlus[misc::TimeSteps][tensor::QInterpolated::size()] = {{}};
  alignas(Alignment) real qInterpolatedMinus[misc::TimeSteps][tensor::QInterpolated::size()] = {{}};
  using QInterpolatedShapeT = real(*)[misc::NumQuantities][misc::NumPaddedPoints];
  auto* qIPlus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedPlus);
  auto* qIMinus = reinterpret_cast<QInterpolatedShapeT>(qInterpolatedMinus);
  // the minus side moves faster and carries more traction, so that no term cancels another
  auto qP = [](std::size_t o, std::size_t q, std::size_t p) {
    return static_cast<real>(1.0 + o + 0.1 * q + 0.01 * p);
  };
  auto qM = [](std::size_t o, std::size_t q, std::size_t p) {
    return static_cast<real>(2.0 + 2.0 * o + 0.3 * q + 0.02 * p);
  };
  for (std::size_t o = 0; o < misc::TimeSteps; ++o) {
    for (std::size_t q = 0; q < misc::NumQuantities; ++q) {
      for (std::size_t p = 0; p < misc::NumPaddedPoints; ++p) {
        qIPlus[o][q][p] = qP(o, q, p);
        qIMinus[o][q][p] = qM(o, q, p);
      }
    }
  }

  ImpedancesAndEta impAndEta{};
  impAndEta.etaS = 7.5e6;
  impAndEta.invZs = 1.0 / 1.5e7;
  impAndEta.invZsNeig = 1.0 / 2.5e7;

  DRGodunovData godunovData{};
  godunovData.doubledSurfaceArea = 3.0;

  real timeWeights[misc::TimeSteps]{};
  for (std::size_t o = 0; o < misc::TimeSteps; ++o) {
    timeWeights[o] = static_cast<real>(0.1 * (o + 1));
  }
  real spaceWeights[misc::NumPaddedPoints]{};
  for (std::size_t p = 0; p < misc::NumPaddedPoints; ++p) {
    spaceWeights[p] = static_cast<real>(0.5 + 0.01 * p);
  }

  // tau* of one node, written out row by row: eta (Y+ tau+ + Y- tau-) where the impedance belongs
  // to the node, the scalar pair of the face everywhere else
  auto interpolated = [&](std::size_t o, std::size_t p, std::size_t row) {
    const double plusTraction[3] = {qP(o, N, p), qP(o, T1, p), qP(o, T2, p)};
    const double minusTraction[3] = {qM(o, N, p), qM(o, T1, p), qM(o, T2, p)};
    if constexpr (NodalImpedance) {
      real plus[9];
      real minus[9];
      real eta[9];
      nodeMatrices(o, p, plus, minus, eta);
      double result = 0.0;
      for (std::size_t k = 0; k < 3; ++k) {
        double admitted = 0.0;
        for (std::size_t c = 0; c < 3; ++c) {
          admitted += plus[3 * c + k] * plusTraction[c] + minus[3 * c + k] * minusTraction[c];
        }
        result += eta[3 * k + row] * admitted;
      }
      return result;
    } else {
      return impAndEta.etaS *
             (impAndEta.invZs * plusTraction[row] + impAndEta.invZsNeig * minusTraction[row]);
    }
  };

  DREnergyOutput energy{};
  FaultStresses<Executor::Host> faultStresses{};
  for (std::size_t o = 0; o < misc::TimeSteps; ++o) {
    fillHost(faultStresses, o);
    common::computeFrictionEnergy(energy,
                                  faultStresses,
                                  qInterpolatedPlus,
                                  qInterpolatedMinus,
                                  impAndEta,
                                  timeWeights[o],
                                  spaceWeights,
                                  godunovData,
                                  o);
  }

  DREnergyOutput deviceEnergy{};
  for (std::uint32_t p = 0; p < misc::NumPaddedPoints; ++p) {
    FaultStresses<Executor::Device> deviceStresses{};
    for (std::uint32_t o = 0; o < misc::TimeSteps; ++o) {
      fillDevice(deviceStresses, o, p);
      common::computeFrictionEnergy<GpuRange>(deviceEnergy,
                                              deviceStresses,
                                              qInterpolatedPlus,
                                              qInterpolatedMinus,
                                              impAndEta,
                                              timeWeights[o],
                                              spaceWeights,
                                              godunovData,
                                              o,
                                              p);
    }
  }

  constexpr double Epsilon = 1e3 * std::numeric_limits<real>::epsilon();
  for (std::size_t p = 0; p < misc::NumPaddedPoints; ++p) {
    double expected = 0.0;
    for (std::size_t o = 0; o < misc::TimeSteps; ++o) {
      const double slipRate2 = qM(o, V, p) - qP(o, V, p);
      const double slipRate3 = qM(o, W, p) - qP(o, W, p);
      const double weight = static_cast<double>(timeWeights[o]) *
                            spaceWeights[p / multisim::NumSimulations] *
                            -godunovData.doubledSurfaceArea;
      expected += weight * (interpolated(o, p, 1) * slipRate2 + interpolated(o, p, 2) * slipRate3);
    }
    CHECK(energy.frictionalEnergy[p] == doctest::Approx(expected).epsilon(Epsilon));
    CHECK(deviceEnergy.frictionalEnergy[p] == doctest::Approx(expected).epsilon(Epsilon));
  }
}

} // namespace seissol::unit_test

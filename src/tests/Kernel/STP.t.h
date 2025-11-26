// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <Equations/Datastructures.h>
#include <Equations/poroelastic/Model/PoroelasticSetup.h>
#include <Model/CommonDatastructures.h>
#include <iomanip>
#include <iostream>

#include <cmath>
#include <limits>
#include <random>
#include <type_traits>

#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"

#include "Equations/poroelastic/Model/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"

namespace seissol::unit_test {

template <typename Cfg>
class SpaceTimePredictorTestFixture {
  public:
  using real = Real<Cfg>;
  constexpr static const double Epsilon = std::numeric_limits<real>::epsilon();
  constexpr static const double Dt = 1.05109e-06;
  real starMatrices0[tensor::star<Cfg>::size(0)];
  real starMatrices1[tensor::star<Cfg>::size(1)];
  real starMatrices2[tensor::star<Cfg>::size(2)];
  real sourceMatrix[tensor::ET<Cfg>::size()];
  real zMatrix[seissol::model::MaterialTT<Cfg>::NumQuantities][tensor::Zinv<Cfg>::size(0)];

  void setStarMatrix(
      const real* at, const real* bt, const real* ct, const double grad[3], real* starMatrix) {
    for (std::size_t idx = 0; idx < seissol::tensor::star<Cfg>::size(0); ++idx) {
      starMatrix[idx] = grad[0] * at[idx];
    }

    for (std::size_t idx = 0; idx < seissol::tensor::star<Cfg>::size(1); ++idx) {
      starMatrix[idx] += grad[1] * bt[idx];
    }

    for (std::size_t idx = 0; idx < seissol::tensor::star<Cfg>::size(2); ++idx) {
      starMatrix[idx] += grad[2] * ct[idx];
    }
  }

  void prepareModel() {
    // prepare Material
    const auto materialVals =
        std::vector<double>{{40.0e9, 2500, 12.0e9, 10.0e9, 0.2, 600.0e-15, 3, 2.5e9, 1040, 0.001}};
    const model::PoroElasticMaterial material(materialVals);

    // prepare Geometry
    std::srand(0);
    std::mt19937 generator(20210109); // Standard mersenne_twister_engine seeded with today's date
    std::uniform_real_distribution<real> distribution(0, 1);
    double x[] = {distribution(generator),
                  distribution(generator),
                  distribution(generator),
                  distribution(generator)};
    double y[] = {distribution(generator),
                  distribution(generator),
                  distribution(generator),
                  distribution(generator)};
    double z[] = {distribution(generator),
                  distribution(generator),
                  distribution(generator),
                  distribution(generator)};
    double gradXi[3];
    double gradEta[3];
    double gradZeta[3];

    transformations::tetrahedronGlobalToReferenceJacobian(x, y, z, gradXi, gradEta, gradZeta);

    // prepare starmatrices
    real atData[tensor::star<Cfg>::size(0)];
    real btData[tensor::star<Cfg>::size(1)];
    real ctData[tensor::star<Cfg>::size(2)];
    auto at = init::star<Cfg>::template view<0>::create(atData);
    auto bt = init::star<Cfg>::template view<0>::create(btData);
    auto ct = init::star<Cfg>::template view<0>::create(ctData);
    model::getTransposedCoefficientMatrix<Cfg>(material, 0, at);
    model::getTransposedCoefficientMatrix<Cfg>(material, 1, bt);
    model::getTransposedCoefficientMatrix<Cfg>(material, 2, ct);
    setStarMatrix(atData, btData, ctData, gradXi, starMatrices0);
    setStarMatrix(atData, btData, ctData, gradEta, starMatrices1);
    setStarMatrix(atData, btData, ctData, gradZeta, starMatrices2);

    // prepare sourceterm
    auto et = init::ET<Cfg>::view::create(sourceMatrix);
    model::getTransposedSourceCoefficientTensor<Cfg>(material, et);

    // prepare Zinv
    auto zinv = init::Zinv<Cfg>::template view<0>::create(zMatrix[0]);
    model::calcZinv<Cfg>(zinv, et, 0, Dt);
    auto zinv1 = init::Zinv<Cfg>::template view<1>::create(zMatrix[1]);
    model::calcZinv<Cfg>(zinv1, et, 1, Dt);
    auto zinv2 = init::Zinv<Cfg>::template view<2>::create(zMatrix[2]);
    model::calcZinv<Cfg>(zinv2, et, 2, Dt);
    auto zinv3 = init::Zinv<Cfg>::template view<3>::create(zMatrix[3]);
    model::calcZinv<Cfg>(zinv3, et, 3, Dt);
    auto zinv4 = init::Zinv<Cfg>::template view<4>::create(zMatrix[4]);
    model::calcZinv<Cfg>(zinv4, et, 4, Dt);
    auto zinv5 = init::Zinv<Cfg>::template view<5>::create(zMatrix[5]);
    model::calcZinv<Cfg>(zinv5, et, 5, Dt);
    auto zinv6 = init::Zinv<Cfg>::template view<6>::create(zMatrix[6]);
    model::calcZinv<Cfg>(zinv6, et, 6, Dt);
    auto zinv7 = init::Zinv<Cfg>::template view<7>::create(zMatrix[7]);
    model::calcZinv<Cfg>(zinv7, et, 7, Dt);
    auto zinv8 = init::Zinv<Cfg>::template view<8>::create(zMatrix[8]);
    model::calcZinv<Cfg>(zinv8, et, 8, Dt);
    auto zinv9 = init::Zinv<Cfg>::template view<9>::create(zMatrix[9]);
    model::calcZinv<Cfg>(zinv9, et, 9, Dt);
    auto zinv10 = init::Zinv<Cfg>::template view<10>::create(zMatrix[10]);
    model::calcZinv<Cfg>(zinv10, et, 10, Dt);
    auto zinv11 = init::Zinv<Cfg>::template view<11>::create(zMatrix[11]);
    model::calcZinv<Cfg>(zinv11, et, 11, Dt);
    auto zinv12 = init::Zinv<Cfg>::template view<12>::create(zMatrix[12]);
    model::calcZinv<Cfg>(zinv12, et, 12, Dt);
  }

  void prepareKernel(seissol::kernel::spaceTimePredictor<Cfg>& krnlPrototype) {
    for (std::size_t n = 0; n < Cfg::ConvergenceOrder; ++n) {
      if (n > 0) {
        for (int d = 0; d < 3; ++d) {
          krnlPrototype.kDivMTSub(d, n) =
              seissol::init::kDivMTSub<Cfg>::Values[seissol::tensor::kDivMTSub<Cfg>::index(d, n)];
        }
      }
      krnlPrototype.selectModes(n) =
          seissol::init::selectModes<Cfg>::Values[seissol::tensor::selectModes<Cfg>::index(n)];
    }
    for (std::size_t k = 0; k < seissol::model::MaterialTT<Cfg>::NumQuantities; k++) {
      krnlPrototype.selectQuantity(k) =
          seissol::init::selectQuantity<Cfg>::Values[seissol::tensor::selectQuantity<Cfg>::index(
              k)];
      krnlPrototype.selectQuantityG(k) =
          init::selectQuantityG<Cfg>::Values[tensor::selectQuantityG<Cfg>::index(k)];
    }
    krnlPrototype.timeInt = seissol::init::timeInt<Cfg>::Values;
    krnlPrototype.wHat = seissol::init::wHat<Cfg>::Values;
  }

  void prepareLHS(seissol::kernel::stpTestLhs<Cfg>& krnlPrototype) {
    krnlPrototype.Z = seissol::init::Z<Cfg>::Values;
    krnlPrototype.deltaLarge = seissol::init::deltaLarge<Cfg>::Values;
    krnlPrototype.deltaSmall = seissol::init::deltaSmall<Cfg>::Values;
  }

  void prepareRHS(seissol::kernel::stpTestRhs<Cfg>& krnlPrototype) {
    for (size_t i = 0; i < 3; i++) {
      krnlPrototype.kDivMT(i) =
          seissol::init::kDivMT<Cfg>::Values[seissol::init::kDivMT<Cfg>::index(i)];
    }
    krnlPrototype.wHat = seissol::init::wHat<Cfg>::Values;
  }

  void prepareQ(real* qData) {
    // scale quantities to make it more realistic
    std::array<real, 13> factor = {{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto q = init::Q<Cfg>::view::create(qData);
    std::srand(1234);
    for (std::size_t qi = 0; qi < q.shape(1); ++qi) {
      for (std::size_t bf = 0; bf < q.shape(0); ++bf) {
        q(bf, qi) = static_cast<real>(std::rand()) / RAND_MAX * factor.at(qi);
      }
    }
  }

  void solveWithKernel(real stp[], const real* qData) {
    real timeIntegrated[seissol::tensor::I<Cfg>::size()];
    alignas(PagesizeStack) real stpRhs[seissol::tensor::spaceTimePredictorRhs<Cfg>::size()];
    std::fill(std::begin(stpRhs), std::end(stpRhs), 0);

    seissol::kernel::spaceTimePredictor<Cfg> krnl;
    prepareKernel(krnl);

    real aValues[seissol::tensor::star<Cfg>::size(0)] = {0};
    real bValues[seissol::tensor::star<Cfg>::size(0)] = {0};
    real cValues[seissol::tensor::star<Cfg>::size(0)] = {0};
    for (size_t i = 0; i < seissol::tensor::star<Cfg>::size(0); i++) {
      aValues[i] = starMatrices0[i] * Dt;
      bValues[i] = starMatrices1[i] * Dt;
      cValues[i] = starMatrices2[i] * Dt;
    }

    krnl.star(0) = aValues;
    krnl.star(1) = bValues;
    krnl.star(2) = cValues;

    for (size_t i = 0; i < seissol::model::MaterialTT<Cfg>::NumQuantities; i++) {
      krnl.Zinv(i) = zMatrix[i];
    }

    auto sourceView = init::ET<Cfg>::view::create(sourceMatrix);
    krnl.Gk = sourceView(10, 6) * Dt;
    krnl.Gl = sourceView(11, 7) * Dt;
    krnl.Gm = sourceView(12, 8) * Dt;

    krnl.Q = qData;
    krnl.I = timeIntegrated;
    krnl.timestep = Dt;
    krnl.spaceTimePredictor = stp;
    krnl.spaceTimePredictorRhs = stpRhs;
    krnl.execute();
  }

  void computeLhs(const real* stp, real* lhs) {
    kernel::stpTestLhs<Cfg> testLhsKrnl;
    prepareLHS(testLhsKrnl);
    testLhsKrnl.ET = sourceMatrix;
    testLhsKrnl.spaceTimePredictor = stp;
    testLhsKrnl.testLhs = lhs;
    testLhsKrnl.minus = -Dt;
    testLhsKrnl.execute();
  };

  void computeRhs(const real* stp, const real* qData, real* rhs) {
    kernel::stpTestRhs<Cfg> testRhsKrnl;
    prepareRHS(testRhsKrnl);
    testRhsKrnl.Q = qData;
    testRhsKrnl.star(0) = starMatrices0;
    testRhsKrnl.star(1) = starMatrices1;
    testRhsKrnl.star(2) = starMatrices2;
    testRhsKrnl.spaceTimePredictor = stp;
    testRhsKrnl.minus = -Dt;
    testRhsKrnl.testRhs = rhs;
    testRhsKrnl.execute();
  };

  SpaceTimePredictorTestFixture() { prepareModel(); };
};

TEST_CASE_TEMPLATE_DEFINE("Solve Space Time Predictor", Cfg, configId1) {
  if constexpr (model::MaterialTT<Cfg>::Type == model::MaterialType::Poroelastic) {
    using real = Real<Cfg>;
    SpaceTimePredictorTestFixture<Cfg> fixture;
    alignas(PagesizeStack) real stp[seissol::tensor::spaceTimePredictor<Cfg>::size()];
    alignas(PagesizeStack) real rhs[seissol::tensor::testLhs<Cfg>::size()];
    alignas(PagesizeStack) real lhs[seissol::tensor::testRhs<Cfg>::size()];
    alignas(PagesizeStack) real qData[seissol::tensor::Q<Cfg>::size()];
    std::fill(std::begin(stp), std::end(stp), 0);
    std::fill(std::begin(rhs), std::end(rhs), 0);
    std::fill(std::begin(lhs), std::end(lhs), 0);
    std::fill(std::begin(qData), std::end(qData), 0);

    fixture.prepareQ(qData);

    fixture.solveWithKernel(stp, qData);

    fixture.computeLhs(stp, lhs);
    fixture.computeRhs(stp, qData, rhs);

    double diffNorm = 0;
    double refNorm = 0;

    auto lhsView = init::testLhs<Cfg>::view::create(lhs);
    auto rhsView = init::testRhs<Cfg>::view::create(rhs);

    for (size_t b = 0; b < tensor::spaceTimePredictor<Cfg>::Shape[0]; b++) {
      for (size_t q = 0; q < tensor::spaceTimePredictor<Cfg>::Shape[1]; q++) {
        for (size_t o = 0; o < tensor::spaceTimePredictor<Cfg>::Shape[2]; o++) {
          const double d = std::abs(lhsView(b, q, o) - rhsView(b, q, o));
          const double a = std::abs(lhsView(b, q, o));
          diffNorm += d * d;
          refNorm += a * a;
        }
      }
    }

    REQUIRE(diffNorm / refNorm < SpaceTimePredictorTestFixture<Cfg>::Epsilon);
  }
}

} // namespace seissol::unit_test

// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <iomanip>
#include <iostream>

#include <cmath>
#include <limits>
#include <random>
#include <type_traits>

#include "Model/Common.h"
#include "Model/PoroelasticSetup.h"
#include "Numerical/Transformation.h"
#include "generated_code/init.h"
#include "generated_code/kernel.h"

#include "Equations/poroelastic/Model/Datastructures.h"
#include "Kernels/Common.h"
#include "generated_code/tensor.h"

namespace seissol::unit_test {

class SpaceTimePredictorTestFixture {
  protected:
  constexpr static auto N =
      seissol::model::MaterialT::NumQuantities * NumBasisFunctions * ConvergenceOrder;
  constexpr static const double Epsilon = std::numeric_limits<real>::epsilon();
  constexpr static const double Dt = 1.05109e-06;
  real starMatrices0[tensor::star::size(0)];
  real starMatrices1[tensor::star::size(1)];
  real starMatrices2[tensor::star::size(2)];
  real sourceMatrix[tensor::ET::size()];
  real zMatrix[seissol::model::MaterialT::NumQuantities][tensor::Zinv::size(0)];

  void setStarMatrix(
      const real* at, const real* bt, const real* ct, const real grad[3], real* starMatrix) {
    for (unsigned idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
      starMatrix[idx] = grad[0] * at[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
      starMatrix[idx] += grad[1] * bt[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
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
    real x[] = {distribution(generator),
                distribution(generator),
                distribution(generator),
                distribution(generator)};
    real y[] = {distribution(generator),
                distribution(generator),
                distribution(generator),
                distribution(generator)};
    real z[] = {distribution(generator),
                distribution(generator),
                distribution(generator),
                distribution(generator)};
    real gradXi[3];
    real gradEta[3];
    real gradZeta[3];

    transformations::tetrahedronGlobalToReferenceJacobian(x, y, z, gradXi, gradEta, gradZeta);

    // prepare starmatrices
    real atData[tensor::star::size(0)];
    real btData[tensor::star::size(1)];
    real ctData[tensor::star::size(2)];
    auto at = init::star::view<0>::create(atData);
    auto bt = init::star::view<0>::create(btData);
    auto ct = init::star::view<0>::create(ctData);
    model::getTransposedCoefficientMatrix(material, 0, at);
    model::getTransposedCoefficientMatrix(material, 1, bt);
    model::getTransposedCoefficientMatrix(material, 2, ct);
    setStarMatrix(atData, btData, ctData, gradXi, starMatrices0);
    setStarMatrix(atData, btData, ctData, gradEta, starMatrices1);
    setStarMatrix(atData, btData, ctData, gradZeta, starMatrices2);

    // prepare sourceterm
    auto et = init::ET::view::create(sourceMatrix);
    model::getTransposedSourceCoefficientTensor(material, et);

    // prepare Zinv
    auto zinv = init::Zinv::view<0>::create(zMatrix[0]);
    model::calcZinv(zinv, et, 0, Dt);
    auto zinv1 = init::Zinv::view<1>::create(zMatrix[1]);
    model::calcZinv(zinv1, et, 1, Dt);
    auto zinv2 = init::Zinv::view<2>::create(zMatrix[2]);
    model::calcZinv(zinv2, et, 2, Dt);
    auto zinv3 = init::Zinv::view<3>::create(zMatrix[3]);
    model::calcZinv(zinv3, et, 3, Dt);
    auto zinv4 = init::Zinv::view<4>::create(zMatrix[4]);
    model::calcZinv(zinv4, et, 4, Dt);
    auto zinv5 = init::Zinv::view<5>::create(zMatrix[5]);
    model::calcZinv(zinv5, et, 5, Dt);
    auto zinv6 = init::Zinv::view<6>::create(zMatrix[6]);
    model::calcZinv(zinv6, et, 6, Dt);
    auto zinv7 = init::Zinv::view<7>::create(zMatrix[7]);
    model::calcZinv(zinv7, et, 7, Dt);
    auto zinv8 = init::Zinv::view<8>::create(zMatrix[8]);
    model::calcZinv(zinv8, et, 8, Dt);
    auto zinv9 = init::Zinv::view<9>::create(zMatrix[9]);
    model::calcZinv(zinv9, et, 9, Dt);
    auto zinv10 = init::Zinv::view<10>::create(zMatrix[10]);
    model::calcZinv(zinv10, et, 10, Dt);
    auto zinv11 = init::Zinv::view<11>::create(zMatrix[11]);
    model::calcZinv(zinv11, et, 11, Dt);
    auto zinv12 = init::Zinv::view<12>::create(zMatrix[12]);
    model::calcZinv(zinv12, et, 12, Dt);
  }

  void prepareKernel(seissol::kernel::spaceTimePredictor& krnlPrototype) {
    for (std::size_t n = 0; n < ConvergenceOrder; ++n) {
      if (n > 0) {
        for (int d = 0; d < 3; ++d) {
          krnlPrototype.kDivMTSub(d, n) =
              seissol::init::kDivMTSub::Values[seissol::tensor::kDivMTSub::index(d, n)];
        }
      }
      krnlPrototype.selectModes(n) =
          seissol::init::selectModes::Values[seissol::tensor::selectModes::index(n)];
    }
    for (std::size_t k = 0; k < seissol::model::MaterialT::NumQuantities; k++) {
      krnlPrototype.selectQuantity(k) =
          seissol::init::selectQuantity::Values[seissol::tensor::selectQuantity::index(k)];
      krnlPrototype.selectQuantityG(k) =
          init::selectQuantityG::Values[tensor::selectQuantityG::index(k)];
    }
    krnlPrototype.timeInt = seissol::init::timeInt::Values;
    krnlPrototype.wHat = seissol::init::wHat::Values;
  }

  void prepareLHS(seissol::kernel::stpTestLhs& krnlPrototype) {
    krnlPrototype.Z = seissol::init::Z::Values;
    krnlPrototype.deltaLarge = seissol::init::deltaLarge::Values;
    krnlPrototype.deltaSmall = seissol::init::deltaSmall::Values;
  }

  void prepareRHS(seissol::kernel::stpTestRhs& krnlPrototype) {
    for (size_t i = 0; i < 3; i++) {
      krnlPrototype.kDivMT(i) = seissol::init::kDivMT::Values[seissol::init::kDivMT::index(i)];
    }
    krnlPrototype.wHat = seissol::init::wHat::Values;
  }

  void prepareQ(real* qData) {
    // scale quantities to make it more realistic
    std::array<real, 13> factor = {{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto q = init::Q::view::create(qData);
    std::srand(1234);
    for (std::size_t qi = 0; qi < seissol::model::MaterialT::NumQuantities; ++qi) {
      for (std::size_t bf = 0; bf < NumBasisFunctions; ++bf) {
        q(bf, qi) = (real)std::rand() / RAND_MAX * factor.at(qi);
      }
    }
  }

  void solveWithKernel(real stp[], const real* qData) {
    real timeIntegrated[seissol::tensor::I::size()];
    alignas(PagesizeStack) real stpRhs[seissol::tensor::spaceTimePredictorRhs::size()];
    std::fill(std::begin(stpRhs), std::end(stpRhs), 0);

    seissol::kernel::spaceTimePredictor krnl;
    prepareKernel(krnl);

    real aValues[seissol::tensor::star::size(0)] = {0};
    real bValues[seissol::tensor::star::size(0)] = {0};
    real cValues[seissol::tensor::star::size(0)] = {0};
    for (size_t i = 0; i < seissol::tensor::star::size(0); i++) {
      aValues[i] = starMatrices0[i] * Dt;
      bValues[i] = starMatrices1[i] * Dt;
      cValues[i] = starMatrices2[i] * Dt;
    }

    krnl.star(0) = aValues;
    krnl.star(1) = bValues;
    krnl.star(2) = cValues;

    for (size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
      krnl.Zinv(i) = zMatrix[i];
    }

    auto sourceView = init::ET::view::create(sourceMatrix);
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
    kernel::stpTestLhs testLhsKrnl;
    prepareLHS(testLhsKrnl);
    testLhsKrnl.ET = sourceMatrix;
    testLhsKrnl.spaceTimePredictor = stp;
    testLhsKrnl.testLhs = lhs;
    testLhsKrnl.minus = -Dt;
    testLhsKrnl.execute();
  };

  void computeRhs(const real* stp, const real* qData, real* rhs) {
    kernel::stpTestRhs testRhsKrnl;
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

  public:
  SpaceTimePredictorTestFixture() { prepareModel(); };
};

TEST_CASE_FIXTURE(SpaceTimePredictorTestFixture, "Solve Space Time Predictor") {
  alignas(PagesizeStack) real stp[seissol::tensor::spaceTimePredictor::size()];
  alignas(PagesizeStack) real rhs[seissol::tensor::testLhs::size()];
  alignas(PagesizeStack) real lhs[seissol::tensor::testRhs::size()];
  alignas(PagesizeStack) real qData[seissol::model::MaterialT::NumQuantities * NumBasisFunctions];
  std::fill(std::begin(stp), std::end(stp), 0);
  std::fill(std::begin(rhs), std::end(rhs), 0);
  std::fill(std::begin(lhs), std::end(lhs), 0);
  std::fill(std::begin(qData), std::end(qData), 0);

  prepareQ(qData);

  solveWithKernel(stp, qData);

  computeLhs(stp, lhs);
  computeRhs(stp, qData, rhs);

  double diffNorm = 0;
  double refNorm = 0;

  auto lhsView = init::testLhs::view::create(lhs);
  auto rhsView = init::testRhs::view::create(rhs);

  for (size_t b = 0; b < tensor::spaceTimePredictor::Shape[0]; b++) {
    for (size_t q = 0; q < tensor::spaceTimePredictor::Shape[1]; q++) {
      for (size_t o = 0; o < tensor::spaceTimePredictor::Shape[2]; o++) {
        const double d = std::abs(lhsView(b, q, o) - rhsView(b, q, o));
        const double a = std::abs(lhsView(b, q, o));
        diffNorm += d * d;
        refNorm += a * a;
      }
    }
  }

  REQUIRE(diffNorm / refNorm < Epsilon);
}

} // namespace seissol::unit_test

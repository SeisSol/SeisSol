// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Real.h"
#include "Config.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"
#include "Equations/poroelastic/Model/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Kernels/Common.h"
#include "Kernels/STP/Setup.h"
#include "Model/Common.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <tuple>
#include <type_traits>

namespace seissol::unit_test {

// The space-time predictor advances the poroelastic configurations; the test runs in the first one.
#define SEISSOL_STP_CONFIG(Cfg) Cfg,
using StpConfig =
    std::tuple_element_t<0, std::tuple<SEISSOL_FOR_EACH_CONFIG_STP(SEISSOL_STP_CONFIG) void>>;
#undef SEISSOL_STP_CONFIG

class SpaceTimePredictorTestFixture {
  protected:
  constexpr static const double Epsilon = std::numeric_limits<Real<StpConfig>>::epsilon();
  constexpr static const double Dt = 1.05109e-06;
  Real<StpConfig> starMatrices0[tensor::star<StpConfig>::size(0)];
  Real<StpConfig> starMatrices1[tensor::star<StpConfig>::size(1)];
  Real<StpConfig> starMatrices2[tensor::star<StpConfig>::size(2)];
  Real<StpConfig> sourceMatrix[tensor::ET<StpConfig>::size()];
  Real<StpConfig> zMatrix[seissol::model::PoroElasticMaterial::NumQuantities]
                         [tensor::Zinv<StpConfig>::size(0)];

  void setStarMatrix(const Real<StpConfig>* at,
                     const Real<StpConfig>* bt,
                     const Real<StpConfig>* ct,
                     const std::array<double, Cell::Dim>& grad,
                     Real<StpConfig>* starMatrix) {
    for (unsigned idx = 0; idx < seissol::tensor::star<StpConfig>::size(0); ++idx) {
      starMatrix[idx] = grad[0] * at[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star<StpConfig>::size(1); ++idx) {
      starMatrix[idx] += grad[1] * bt[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star<StpConfig>::size(2); ++idx) {
      starMatrix[idx] += grad[2] * ct[idx];
    }
  }

  void prepareModel() {
    // prepare Material
    const auto materialVals =
        std::vector<double>{{40.0e9, 2500, 12.0e9, 10.0e9, 0.2, 600.0e-15, 3, 2.5e9, 1040, 0.001}};
    const model::PoroElasticMaterial material(materialVals);

    // NOLINTNEXTLINE (-cert-dcl59-cpp)
    std::mt19937 generator(20210109); // Standard mersenne_twister_engine seeded with today's date
    std::uniform_real_distribution<Real<StpConfig>> distribution(0, 1);
    std::array<CoordinateT, Cell::NumVertices> vertices{};
    for (auto& vertex : vertices) {
      vertex =
          CoordinateT{distribution(generator), distribution(generator), distribution(generator)};
    }

    const auto grad = seissol::geometry::AffineTransform(vertices).refToSpaceJacobianInverse(
        seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data()));

    std::array<double, Cell::Dim> gradXi{};
    std::array<double, Cell::Dim> gradEta{};
    std::array<double, Cell::Dim> gradZeta{};
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      gradXi[i] = grad(0, i);
      gradEta[i] = grad(1, i);
      gradZeta[i] = grad(2, i);
    }

    // prepare starmatrices
    Real<StpConfig> atData[tensor::star<StpConfig>::size(0)];
    Real<StpConfig> btData[tensor::star<StpConfig>::size(1)];
    Real<StpConfig> ctData[tensor::star<StpConfig>::size(2)];
    auto at = init::star<StpConfig>::view<0>::create(atData);
    auto bt = init::star<StpConfig>::view<0>::create(btData);
    auto ct = init::star<StpConfig>::view<0>::create(ctData);
    model::getTransposedCoefficientMatrix<StpConfig>(material, 0, at);
    model::getTransposedCoefficientMatrix<StpConfig>(material, 1, bt);
    model::getTransposedCoefficientMatrix<StpConfig>(material, 2, ct);
    setStarMatrix(atData, btData, ctData, gradXi, starMatrices0);
    setStarMatrix(atData, btData, ctData, gradEta, starMatrices1);
    setStarMatrix(atData, btData, ctData, gradZeta, starMatrices2);

    // prepare sourceterm
    auto et = init::ET<StpConfig>::view::create(sourceMatrix);
    model::getTransposedSourceCoefficientTensor<StpConfig>(material, et);

    // prepare Zinv
    auto zinv0 = init::Zinv<StpConfig>::view<0>::create(zMatrix[0]);
    model::calcZinv<StpConfig>(zinv0, et, 0, model::isStiffRow<model::PoroElasticMaterial>(0), Dt);
    auto zinv1 = init::Zinv<StpConfig>::view<1>::create(zMatrix[1]);
    model::calcZinv<StpConfig>(zinv1, et, 1, model::isStiffRow<model::PoroElasticMaterial>(1), Dt);
    auto zinv2 = init::Zinv<StpConfig>::view<2>::create(zMatrix[2]);
    model::calcZinv<StpConfig>(zinv2, et, 2, model::isStiffRow<model::PoroElasticMaterial>(2), Dt);
    auto zinv3 = init::Zinv<StpConfig>::view<3>::create(zMatrix[3]);
    model::calcZinv<StpConfig>(zinv3, et, 3, model::isStiffRow<model::PoroElasticMaterial>(3), Dt);
    auto zinv4 = init::Zinv<StpConfig>::view<4>::create(zMatrix[4]);
    model::calcZinv<StpConfig>(zinv4, et, 4, model::isStiffRow<model::PoroElasticMaterial>(4), Dt);
    auto zinv5 = init::Zinv<StpConfig>::view<5>::create(zMatrix[5]);
    model::calcZinv<StpConfig>(zinv5, et, 5, model::isStiffRow<model::PoroElasticMaterial>(5), Dt);
    auto zinv6 = init::Zinv<StpConfig>::view<6>::create(zMatrix[6]);
    model::calcZinv<StpConfig>(zinv6, et, 6, model::isStiffRow<model::PoroElasticMaterial>(6), Dt);
    auto zinv7 = init::Zinv<StpConfig>::view<7>::create(zMatrix[7]);
    model::calcZinv<StpConfig>(zinv7, et, 7, model::isStiffRow<model::PoroElasticMaterial>(7), Dt);
    auto zinv8 = init::Zinv<StpConfig>::view<8>::create(zMatrix[8]);
    model::calcZinv<StpConfig>(zinv8, et, 8, model::isStiffRow<model::PoroElasticMaterial>(8), Dt);
    auto zinv9 = init::Zinv<StpConfig>::view<9>::create(zMatrix[9]);
    model::calcZinv<StpConfig>(zinv9, et, 9, model::isStiffRow<model::PoroElasticMaterial>(9), Dt);
    auto zinv10 = init::Zinv<StpConfig>::view<10>::create(zMatrix[10]);
    model::calcZinv<StpConfig>(
        zinv10, et, 10, model::isStiffRow<model::PoroElasticMaterial>(10), Dt);
    auto zinv11 = init::Zinv<StpConfig>::view<11>::create(zMatrix[11]);
    model::calcZinv<StpConfig>(
        zinv11, et, 11, model::isStiffRow<model::PoroElasticMaterial>(11), Dt);
    auto zinv12 = init::Zinv<StpConfig>::view<12>::create(zMatrix[12]);
    model::calcZinv<StpConfig>(
        zinv12, et, 12, model::isStiffRow<model::PoroElasticMaterial>(12), Dt);
  }

  void prepareKernel(seissol::kernel::spaceTimePredictor<StpConfig>& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool<StpConfig>::host());
  }

  void prepareLHS(seissol::kernel::stpTestLhs<StpConfig>& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool<StpConfig>::host());
  }

  void prepareRHS(seissol::kernel::stpTestRhs<StpConfig>& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool<StpConfig>::host());
  }

  void prepareQ(Real<StpConfig>* qData) {
    // scale quantities to make it more realistic
    std::array<Real<StpConfig>, 13> factor = {
        {1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto q = init::Q<StpConfig>::view::create(qData);

    // NOLINTNEXTLINE (-cert-dcl59-cpp)
    std::mt19937 rnggen(1234);
    std::uniform_real_distribution<> rngdist(0.0, 1.0);

    for (std::size_t qi = 0; qi < q.shape(1); ++qi) {
      for (std::size_t bf = 0; bf < q.shape(0); ++bf) {
        q(bf, qi) = rngdist(rnggen) * factor.at(qi);
      }
    }
  }

  void solveWithKernel(Real<StpConfig> stp[], const Real<StpConfig>* qData) {
    Real<StpConfig> timeIntegrated[seissol::tensor::I<StpConfig>::size()];

    seissol::kernel::spaceTimePredictor<StpConfig> krnl;
    prepareKernel(krnl);

    Real<StpConfig> aValues[seissol::tensor::star<StpConfig>::size(0)] = {0};
    Real<StpConfig> bValues[seissol::tensor::star<StpConfig>::size(0)] = {0};
    Real<StpConfig> cValues[seissol::tensor::star<StpConfig>::size(0)] = {0};

    // Scaled by Dt, as Spacetime::executeSTP does. The minus sign of the flux term is not
    // applied here: kDivMT carries it, negated at code generation (negateFamily in
    // codegen/kernels/aderdg/aderdg.py).
    for (size_t i = 0; i < seissol::tensor::star<StpConfig>::size(0); i++) {
      aValues[i] = starMatrices0[i] * Dt;
      bValues[i] = starMatrices1[i] * Dt;
      cValues[i] = starMatrices2[i] * Dt;
    }

    krnl.star(0) = aValues;
    krnl.star(1) = bValues;
    krnl.star(2) = cValues;

    for (size_t i = 0; i < seissol::model::PoroElasticMaterial::NumQuantities; i++) {
      krnl.Zinv(i) = zMatrix[i];
    }

    auto sourceView = init::ET<StpConfig>::view::create(sourceMatrix);
    for (std::size_t i = 0; i < model::PoroElasticMaterial::StiffSourceRows.size(); ++i) {
      const auto& row = model::PoroElasticMaterial::StiffSourceRows[i];
      krnl.G(i) = sourceView(row.quantity, row.target) * Dt;
    }

    krnl.Q = qData;
    krnl.I = timeIntegrated;
    krnl.timestep = Dt;
    krnl.spaceTimePredictor = stp;
    krnl.execute();
  }

  void computeLhs(const Real<StpConfig>* stp, Real<StpConfig>* lhs) {
    kernel::stpTestLhs<StpConfig> testLhsKrnl;
    prepareLHS(testLhsKrnl);
    testLhsKrnl.ET = sourceMatrix;
    testLhsKrnl.spaceTimePredictor = stp;
    testLhsKrnl.testLhs = lhs;
    testLhsKrnl.minus = -Dt;
    testLhsKrnl.execute();
  };

  void computeRhs(const Real<StpConfig>* stp, const Real<StpConfig>* qData, Real<StpConfig>* rhs) {
    kernel::stpTestRhs<StpConfig> testRhsKrnl;
    prepareRHS(testRhsKrnl);
    testRhsKrnl.Q = qData;
    testRhsKrnl.star(0) = starMatrices0;
    testRhsKrnl.star(1) = starMatrices1;
    testRhsKrnl.star(2) = starMatrices2;
    testRhsKrnl.spaceTimePredictor = stp;
    // The flux term is -Dt * star * K^T; kDivMT already is -K^T, so its factor here is +Dt.
    testRhsKrnl.minus = Dt;
    testRhsKrnl.testRhs = rhs;
    testRhsKrnl.execute();
  };

  public:
  SpaceTimePredictorTestFixture() { prepareModel(); };
};

TEST_CASE_FIXTURE(SpaceTimePredictorTestFixture,
                  "Solve Space Time Predictor" * doctest::test_suite("kernel")) {
  alignas(PagesizeStack) Real<StpConfig>
      stp[seissol::tensor::spaceTimePredictor<StpConfig>::size()];
  alignas(PagesizeStack) Real<StpConfig> rhs[seissol::tensor::testLhs<StpConfig>::size()];
  alignas(PagesizeStack) Real<StpConfig> lhs[seissol::tensor::testRhs<StpConfig>::size()];
  alignas(PagesizeStack) Real<StpConfig> qData[seissol::tensor::Q<StpConfig>::size()];
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

  auto lhsView = init::testLhs<StpConfig>::view::create(lhs);
  auto rhsView = init::testRhs<StpConfig>::view::create(rhs);

  for (size_t b = 0; b < tensor::spaceTimePredictor<StpConfig>::Shape[0]; b++) {
    for (size_t q = 0; q < tensor::spaceTimePredictor<StpConfig>::Shape[1]; q++) {
      for (size_t o = 0; o < tensor::spaceTimePredictor<StpConfig>::Shape[2]; o++) {
        const double d = std::abs(lhsView(b, q, o) - rhsView(b, q, o));
        const double a = std::abs(lhsView(b, q, o));
        diffNorm += d * d;
        refNorm += a * a;
      }
    }
  }

  CHECK(diffNorm / refNorm < Epsilon);
}

} // namespace seissol::unit_test

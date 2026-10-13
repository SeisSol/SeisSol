// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Real.h"
#include "Common/Typedefs.h"
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
#include "TestConfigs.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>

namespace seissol::unit_test {

// The space-time predictor advances the poroelastic configurations; the test runs for each of them.
template <typename Cfg>
class SpaceTimePredictorTestFixture {
  protected:
  constexpr static const double Epsilon = std::numeric_limits<Real<Cfg>>::epsilon();
  constexpr static const double Dt = 1.05109e-06;
  Real<Cfg> starMatrices0_[tensor::star<Cfg>::size(0)];
  Real<Cfg> starMatrices1_[tensor::star<Cfg>::size(1)];
  Real<Cfg> starMatrices2_[tensor::star<Cfg>::size(2)];
  Real<Cfg> sourceMatrix_[tensor::ET<Cfg>::size()];
  Real<Cfg> zMatrix_[seissol::model::PoroElasticMaterial::NumQuantities]
                    [tensor::Zinv<Cfg>::size(0)];

  void setStarMatrix(const Real<Cfg>* at,
                     const Real<Cfg>* bt,
                     const Real<Cfg>* ct,
                     const std::array<double, Cell::Dim>& grad,
                     Real<Cfg>* starMatrix) {
    for (unsigned idx = 0; idx < seissol::tensor::star<Cfg>::size(0); ++idx) {
      starMatrix[idx] = grad[0] * at[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star<Cfg>::size(1); ++idx) {
      starMatrix[idx] += grad[1] * bt[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star<Cfg>::size(2); ++idx) {
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
    std::uniform_real_distribution<Real<Cfg>> distribution(0, 1);
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
    Real<Cfg> atData[tensor::star<Cfg>::size(0)];
    Real<Cfg> btData[tensor::star<Cfg>::size(1)];
    Real<Cfg> ctData[tensor::star<Cfg>::size(2)];
    auto at = init::star<Cfg>::template view<0>::create(atData);
    auto bt = init::star<Cfg>::template view<0>::create(btData);
    auto ct = init::star<Cfg>::template view<0>::create(ctData);
    model::getTransposedCoefficientMatrix<Cfg>(material, 0, at);
    model::getTransposedCoefficientMatrix<Cfg>(material, 1, bt);
    model::getTransposedCoefficientMatrix<Cfg>(material, 2, ct);
    setStarMatrix(atData, btData, ctData, gradXi, starMatrices0_);
    setStarMatrix(atData, btData, ctData, gradEta, starMatrices1_);
    setStarMatrix(atData, btData, ctData, gradZeta, starMatrices2_);

    // prepare sourceterm
    auto et = init::ET<Cfg>::view::create(sourceMatrix_);
    model::getTransposedSourceCoefficientTensor<Cfg>(material, et);

    // prepare Zinv
    auto zinv0 = init::Zinv<Cfg>::template view<0>::create(zMatrix_[0]);
    model::calcZinv<Cfg>(zinv0, et, 0, model::isStiffRow<model::PoroElasticMaterial>(0), Dt);
    auto zinv1 = init::Zinv<Cfg>::template view<1>::create(zMatrix_[1]);
    model::calcZinv<Cfg>(zinv1, et, 1, model::isStiffRow<model::PoroElasticMaterial>(1), Dt);
    auto zinv2 = init::Zinv<Cfg>::template view<2>::create(zMatrix_[2]);
    model::calcZinv<Cfg>(zinv2, et, 2, model::isStiffRow<model::PoroElasticMaterial>(2), Dt);
    auto zinv3 = init::Zinv<Cfg>::template view<3>::create(zMatrix_[3]);
    model::calcZinv<Cfg>(zinv3, et, 3, model::isStiffRow<model::PoroElasticMaterial>(3), Dt);
    auto zinv4 = init::Zinv<Cfg>::template view<4>::create(zMatrix_[4]);
    model::calcZinv<Cfg>(zinv4, et, 4, model::isStiffRow<model::PoroElasticMaterial>(4), Dt);
    auto zinv5 = init::Zinv<Cfg>::template view<5>::create(zMatrix_[5]);
    model::calcZinv<Cfg>(zinv5, et, 5, model::isStiffRow<model::PoroElasticMaterial>(5), Dt);
    auto zinv6 = init::Zinv<Cfg>::template view<6>::create(zMatrix_[6]);
    model::calcZinv<Cfg>(zinv6, et, 6, model::isStiffRow<model::PoroElasticMaterial>(6), Dt);
    auto zinv7 = init::Zinv<Cfg>::template view<7>::create(zMatrix_[7]);
    model::calcZinv<Cfg>(zinv7, et, 7, model::isStiffRow<model::PoroElasticMaterial>(7), Dt);
    auto zinv8 = init::Zinv<Cfg>::template view<8>::create(zMatrix_[8]);
    model::calcZinv<Cfg>(zinv8, et, 8, model::isStiffRow<model::PoroElasticMaterial>(8), Dt);
    auto zinv9 = init::Zinv<Cfg>::template view<9>::create(zMatrix_[9]);
    model::calcZinv<Cfg>(zinv9, et, 9, model::isStiffRow<model::PoroElasticMaterial>(9), Dt);
    auto zinv10 = init::Zinv<Cfg>::template view<10>::create(zMatrix_[10]);
    model::calcZinv<Cfg>(zinv10, et, 10, model::isStiffRow<model::PoroElasticMaterial>(10), Dt);
    auto zinv11 = init::Zinv<Cfg>::template view<11>::create(zMatrix_[11]);
    model::calcZinv<Cfg>(zinv11, et, 11, model::isStiffRow<model::PoroElasticMaterial>(11), Dt);
    auto zinv12 = init::Zinv<Cfg>::template view<12>::create(zMatrix_[12]);
    model::calcZinv<Cfg>(zinv12, et, 12, model::isStiffRow<model::PoroElasticMaterial>(12), Dt);
  }

  void prepareKernel(seissol::kernel::spaceTimePredictor<Cfg>& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool<Cfg>::host());
  }

  void prepareLHS(seissol::kernel::stpTestLhs<Cfg>& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool<Cfg>::host());
  }

  void prepareRHS(seissol::kernel::stpTestRhs<Cfg>& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool<Cfg>::host());
  }

  void prepareQ(Real<Cfg>* qData) {
    // scale quantities to make it more realistic
    std::array<Real<Cfg>, 13> factor = {{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto q = init::Q<Cfg>::view::create(qData);

    // NOLINTNEXTLINE (-cert-dcl59-cpp)
    std::mt19937 rnggen(1234);
    std::uniform_real_distribution<> rngdist(0.0, 1.0);

    for (std::size_t qi = 0; qi < q.shape(1); ++qi) {
      for (std::size_t bf = 0; bf < q.shape(0); ++bf) {
        q(bf, qi) = rngdist(rnggen) * factor.at(qi);
      }
    }
  }

  void solveWithKernel(Real<Cfg> stp[], const Real<Cfg>* qData) {
    Real<Cfg> timeIntegrated[seissol::tensor::I<Cfg>::size()];

    seissol::kernel::spaceTimePredictor<Cfg> krnl;
    prepareKernel(krnl);

    Real<Cfg> aValues[seissol::tensor::star<Cfg>::size(0)] = {0};
    Real<Cfg> bValues[seissol::tensor::star<Cfg>::size(0)] = {0};
    Real<Cfg> cValues[seissol::tensor::star<Cfg>::size(0)] = {0};

    // Scaled by Dt, as Spacetime::executeSTP does. The minus sign of the flux term is not
    // applied here: kDivMT carries it, negated at code generation (negateFamily in
    // codegen/kernels/aderdg/aderdg.py).
    for (size_t i = 0; i < seissol::tensor::star<Cfg>::size(0); i++) {
      aValues[i] = starMatrices0_[i] * Dt;
      bValues[i] = starMatrices1_[i] * Dt;
      cValues[i] = starMatrices2_[i] * Dt;
    }

    krnl.star(0) = aValues;
    krnl.star(1) = bValues;
    krnl.star(2) = cValues;

    for (size_t i = 0; i < seissol::model::PoroElasticMaterial::NumQuantities; i++) {
      krnl.Zinv(i) = zMatrix_[i];
    }

    auto sourceView = init::ET<Cfg>::view::create(sourceMatrix_);
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

  void computeLhs(const Real<Cfg>* stp, Real<Cfg>* lhs) {
    kernel::stpTestLhs<Cfg> testLhsKrnl;
    prepareLHS(testLhsKrnl);
    testLhsKrnl.ET = sourceMatrix_;
    testLhsKrnl.spaceTimePredictor = stp;
    testLhsKrnl.testLhs = lhs;
    testLhsKrnl.minus = -Dt;
    testLhsKrnl.execute();
  };

  void computeRhs(const Real<Cfg>* stp, const Real<Cfg>* qData, Real<Cfg>* rhs) {
    kernel::stpTestRhs<Cfg> testRhsKrnl;
    prepareRHS(testRhsKrnl);
    testRhsKrnl.Q = qData;
    testRhsKrnl.star(0) = starMatrices0_;
    testRhsKrnl.star(1) = starMatrices1_;
    testRhsKrnl.star(2) = starMatrices2_;
    testRhsKrnl.spaceTimePredictor = stp;
    // The flux term is -Dt * star * K^T; kDivMT already is -K^T, so its factor here is +Dt.
    testRhsKrnl.minus = Dt;
    testRhsKrnl.testRhs = rhs;
    testRhsKrnl.execute();
  };

  public:
  SpaceTimePredictorTestFixture() { prepareModel(); };

  void solveAndCompare() {
    alignas(PagesizeStack) Real<Cfg> stp[seissol::tensor::spaceTimePredictor<Cfg>::size()];
    alignas(PagesizeStack) Real<Cfg> rhs[seissol::tensor::testLhs<Cfg>::size()];
    alignas(PagesizeStack) Real<Cfg> lhs[seissol::tensor::testRhs<Cfg>::size()];
    alignas(PagesizeStack) Real<Cfg> qData[seissol::tensor::Q<Cfg>::size()];
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

    CHECK(diffNorm / refNorm < Epsilon);
  }
};

TEST_CASE_TEMPLATE_DEFINE("Solve Space Time Predictor" * doctest::test_suite("kernel"),
                          Cfg,
                          SolveSpaceTimePredictor) {
  SpaceTimePredictorTestFixture<Cfg> fixture;
  fixture.solveAndCompare();
}

TEST_CASE_TEMPLATE_APPLY(SolveSpaceTimePredictor, ConfigsOfSolver<SolverType::STP>);

} // namespace seissol::unit_test

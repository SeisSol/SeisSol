#include <iostream>
#include <iomanip>

#include <limits>
#include <cmath>
#include <type_traits>
#include <random>

#include "generated_code/kernel.h"
#include "generated_code/init.h"
#include <Numerical_aux/Transformation.h>
#include <Model/common.hpp>
#include <Model/PoroelasticSetup.h>

#include "Kernels/common.hpp"
#include "Equations/poroelastic/Model/datastructures.hpp"
#include "generated_code/tensor.h"

namespace seissol::unit_test {

class SpaceTimePredictorTestFixture {
  protected:
  const int N = NUMBER_OF_QUANTITIES * NUMBER_OF_BASIS_FUNCTIONS * CONVERGENCE_ORDER;
  constexpr static double const epsilon = std::numeric_limits<real>::epsilon();
  constexpr static double const dt = 1.05109e-06;
  real starMatrices0[tensor::star::size(0)];
  real starMatrices1[tensor::star::size(1)];
  real starMatrices2[tensor::star::size(2)];
  real sourceMatrix[tensor::ET::size()];
  real zMatrix[NUMBER_OF_QUANTITIES][tensor::Zinv::size(0)];

  void setStarMatrix(real* i_AT, real* i_BT, real* i_CT, real i_grad[3], real* o_starMatrix) {
    for (unsigned idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
      o_starMatrix[idx] = i_grad[0] * i_AT[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
      o_starMatrix[idx] += i_grad[1] * i_BT[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
      o_starMatrix[idx] += i_grad[2] * i_CT[idx];
    }
  }

  void prepareModel() {
    // prepare Material
    std::array<double, 10> materialVals = {
        {40.0e9, 2500, 12.0e9, 10.0e9, 0.2, 600.0e-15, 3, 2.5e9, 1040, 0.001}};
    model::PoroElasticMaterial material(materialVals.data(), 10);

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
    real ATData[tensor::star::size(0)];
    real BTData[tensor::star::size(1)];
    real CTData[tensor::star::size(2)];
    auto AT = init::star::view<0>::create(ATData);
    auto BT = init::star::view<0>::create(BTData);
    auto CT = init::star::view<0>::create(CTData);
    model::getTransposedCoefficientMatrix(material, 0, AT);
    model::getTransposedCoefficientMatrix(material, 1, BT);
    model::getTransposedCoefficientMatrix(material, 2, CT);
    setStarMatrix(ATData, BTData, CTData, gradXi, starMatrices0);
    setStarMatrix(ATData, BTData, CTData, gradEta, starMatrices1);
    setStarMatrix(ATData, BTData, CTData, gradZeta, starMatrices2);

    // prepare sourceterm
    auto ET = init::ET::view::create(sourceMatrix);
    model::getTransposedSourceCoefficientTensor(material, ET);

    // prepare Zinv
    auto Zinv = init::Zinv::view<0>::create(zMatrix[0]);
    model::calcZinv(Zinv, ET, 0, dt);
    auto Zinv1 = init::Zinv::view<1>::create(zMatrix[1]);
    model::calcZinv(Zinv1, ET, 1, dt);
    auto Zinv2 = init::Zinv::view<2>::create(zMatrix[2]);
    model::calcZinv(Zinv2, ET, 2, dt);
    auto Zinv3 = init::Zinv::view<3>::create(zMatrix[3]);
    model::calcZinv(Zinv3, ET, 3, dt);
    auto Zinv4 = init::Zinv::view<4>::create(zMatrix[4]);
    model::calcZinv(Zinv4, ET, 4, dt);
    auto Zinv5 = init::Zinv::view<5>::create(zMatrix[5]);
    model::calcZinv(Zinv5, ET, 5, dt);
    auto Zinv6 = init::Zinv::view<6>::create(zMatrix[6]);
    model::calcZinv(Zinv6, ET, 6, dt);
    auto Zinv7 = init::Zinv::view<7>::create(zMatrix[7]);
    model::calcZinv(Zinv7, ET, 7, dt);
    auto Zinv8 = init::Zinv::view<8>::create(zMatrix[8]);
    model::calcZinv(Zinv8, ET, 8, dt);
    auto Zinv9 = init::Zinv::view<9>::create(zMatrix[9]);
    model::calcZinv(Zinv9, ET, 9, dt);
    auto Zinv10 = init::Zinv::view<10>::create(zMatrix[10]);
    model::calcZinv(Zinv10, ET, 10, dt);
    auto Zinv11 = init::Zinv::view<11>::create(zMatrix[11]);
    model::calcZinv(Zinv11, ET, 11, dt);
    auto Zinv12 = init::Zinv::view<12>::create(zMatrix[12]);
    model::calcZinv(Zinv12, ET, 12, dt);
  }

  void prepareKernel(seissol::kernel::spaceTimePredictor& m_krnlPrototype) {
    for (int n = 0; n < CONVERGENCE_ORDER; ++n) {
      if (n > 0) {
        for (int d = 0; d < 3; ++d) {
          m_krnlPrototype.kDivMTSub(d, n) =
              seissol::init::kDivMTSub::Values[seissol::tensor::kDivMTSub::index(d, n)];
        }
      }
      m_krnlPrototype.selectModes(n) =
          seissol::init::selectModes::Values[seissol::tensor::selectModes::index(n)];
    }
    for (int k = 0; k < NUMBER_OF_QUANTITIES; k++) {
      m_krnlPrototype.selectQuantity(k) =
          seissol::init::selectQuantity::Values[seissol::tensor::selectQuantity::index(k)];
      m_krnlPrototype.selectQuantityG(k) =
          init::selectQuantityG::Values[tensor::selectQuantityG::index(k)];
    }
    m_krnlPrototype.timeInt = seissol::init::timeInt::Values;
    m_krnlPrototype.wHat = seissol::init::wHat::Values;
  }

  void prepareLHS(seissol::kernel::stpTestLhs& m_krnlPrototype) {
    m_krnlPrototype.Z = seissol::init::Z::Values;
    m_krnlPrototype.deltaLarge = seissol::init::deltaLarge::Values;
    m_krnlPrototype.deltaSmall = seissol::init::deltaSmall::Values;
  }

  void prepareRHS(seissol::kernel::stpTestRhs& m_krnlPrototype) {
    for (size_t i = 0; i < 3; i++) {
      m_krnlPrototype.kDivMT(i) = seissol::init::kDivMT::Values[seissol::init::kDivMT::index(i)];
    }
    m_krnlPrototype.wHat = seissol::init::wHat::Values;
  }

  void prepareQ(real* QData) {
    // scale quantities to make it more realistic
    std::array<real, 13> factor = {{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto Q = init::Q::view::create(QData);
    std::srand(1234);
    for (int q = 0; q < NUMBER_OF_QUANTITIES; q++) {
      for (int bf = 0; bf < NUMBER_OF_BASIS_FUNCTIONS; bf++) {
        Q(bf, q) = (real)std::rand() / RAND_MAX * factor.at(q);
      }
    }
  }

  void solveWithKernel(real stp[], real* QData) {
    real o_timeIntegrated[seissol::tensor::I::size()];
    alignas(PAGESIZE_STACK) real stpRhs[seissol::tensor::spaceTimePredictorRhs::size()];
    std::fill(std::begin(stpRhs), std::end(stpRhs), 0);

    seissol::kernel::spaceTimePredictor krnl;
    prepareKernel(krnl);

    real A_values[seissol::tensor::star::size(0)] = {0};
    real B_values[seissol::tensor::star::size(0)] = {0};
    real C_values[seissol::tensor::star::size(0)] = {0};
    for (size_t i = 0; i < seissol::tensor::star::size(0); i++) {
      A_values[i] = starMatrices0[i] * dt;
      B_values[i] = starMatrices1[i] * dt;
      C_values[i] = starMatrices2[i] * dt;
    }

    krnl.star(0) = A_values;
    krnl.star(1) = B_values;
    krnl.star(2) = C_values;

    for (size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
      krnl.Zinv(i) = zMatrix[i];
    }

    auto sourceView = init::ET::view::create(sourceMatrix);
    krnl.Gk = sourceView(10, 6) * dt;
    krnl.Gl = sourceView(11, 7) * dt;
    krnl.Gm = sourceView(12, 8) * dt;

    krnl.Q = QData;
    krnl.I = o_timeIntegrated;
    krnl.timestep = dt;
    krnl.spaceTimePredictor = stp;
    krnl.spaceTimePredictorRhs = stpRhs;
    krnl.execute();
  }

  void computeLhs(real* stp, real* lhs) {
    kernel::stpTestLhs testLhsKrnl;
    prepareLHS(testLhsKrnl);
    testLhsKrnl.ET = sourceMatrix;
    testLhsKrnl.spaceTimePredictor = stp;
    testLhsKrnl.testLhs = lhs;
    testLhsKrnl.minus = -dt;
    testLhsKrnl.execute();
  };

  void computeRhs(real* stp, real* QData, real* rhs) {
    kernel::stpTestRhs testRhsKrnl;
    prepareRHS(testRhsKrnl);
    testRhsKrnl.Q = QData;
    testRhsKrnl.star(0) = starMatrices0;
    testRhsKrnl.star(1) = starMatrices1;
    testRhsKrnl.star(2) = starMatrices2;
    testRhsKrnl.spaceTimePredictor = stp;
    testRhsKrnl.minus = -dt;
    testRhsKrnl.testRhs = rhs;
    testRhsKrnl.execute();
  };

  public:
  SpaceTimePredictorTestFixture() { prepareModel(); };
};

TEST_CASE_FIXTURE(SpaceTimePredictorTestFixture, "Solve Space Time Predictor") {
  alignas(PAGESIZE_STACK) real stp[seissol::tensor::spaceTimePredictor::size()];
  alignas(PAGESIZE_STACK) real rhs[seissol::tensor::testLhs::size()];
  alignas(PAGESIZE_STACK) real lhs[seissol::tensor::testRhs::size()];
  real QData[NUMBER_OF_QUANTITIES * NUMBER_OF_BASIS_FUNCTIONS];
  std::fill(std::begin(stp), std::end(stp), 0);
  std::fill(std::begin(rhs), std::end(rhs), 0);
  std::fill(std::begin(lhs), std::end(lhs), 0);
  std::fill(std::begin(QData), std::end(QData), 0);

  prepareQ(QData);

  solveWithKernel(stp, QData);

  computeLhs(stp, lhs);
  computeRhs(stp, QData, rhs);

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

  REQUIRE(diffNorm / refNorm < epsilon);
}

} // namespace seissol::unit_test
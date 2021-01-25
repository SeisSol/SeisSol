#include <iostream>
#include <iomanip>

#include <cxxtest/TestSuite.h>
#include <limits>
#include <math.h>
#include <type_traits>

#include "generated_code/kernel.h"
#include "generated_code/init.h"
#include "data.h"

#include "Equations/poroelastic/Model/datastructures.hpp"

namespace seissol {
  namespace unit_test {
    class SpaceTimeTestSuite;
  }
}

class seissol::unit_test::SpaceTimeTestSuite : public CxxTest::TestSuite
{
  private:
  const int N = NUMBER_OF_QUANTITIES*NUMBER_OF_BASIS_FUNCTIONS*CONVERGENCE_ORDER;
  double epsilon = std::numeric_limits<real>::epsilon();

  void setStarMatrix( real* i_AT,
      real* i_BT,
      real* i_CT,
      real  i_grad[3],
      real* o_starMatrix )
  {
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
  
  void prepareModel(real* starMatrices0, real*starMatrices1, real* starMatrices2, real* sourceMatrix, real zMatrix[NUMBER_OF_QUANTITIES][CONVERGENCE_ORDER*CONVERGENCE_ORDER]) {
    //prepare Material
    std::array<double, 10> materialVals = {{40.0e9, 2500, 12.0e9, 10.0e9, 0.2, 600.0e-15, 3, 2.5e9, 1040, 0.001}};
    model::PoroElasticMaterial material(materialVals.data(), 10);

    //prepare Geometry
    std::srand(0);
    real x[] = {(real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX};
    real y[] = {(real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX};
    real z[] = {(real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX, (real)std::rand()/RAND_MAX};
    real gradXi[3];
    real gradEta[3];
    real gradZeta[3];

    transformations::tetrahedronGlobalToReferenceJacobian( x, y, z, gradXi, gradEta, gradZeta );

    //prepare starmatrices
    real ATData[tensor::star::size(0)];
    real BTData[tensor::star::size(1)];
    real CTData[tensor::star::size(2)];
    auto AT = init::star::view<0>::create(ATData);
    auto BT = init::star::view<0>::create(BTData);
    auto CT = init::star::view<0>::create(CTData);
    model::getTransposedCoefficientMatrix(material, 0, AT );
    model::getTransposedCoefficientMatrix(material, 1, BT );
    model::getTransposedCoefficientMatrix(material, 2, CT );
    setStarMatrix(ATData, BTData, CTData, gradXi,   starMatrices0);
    setStarMatrix(ATData, BTData, CTData, gradEta,  starMatrices1);
    setStarMatrix(ATData, BTData, CTData, gradZeta, starMatrices2);

    //prepare sourceterm
    auto ET = init::ET::view::create(sourceMatrix);
    model::getTransposedSourceCoefficientTensor(material, ET);

    //prepare Zinv
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

  void prepareKernel(seissol::kernel::stp& m_krnlPrototype) { 
    for (int n = 0; n < CONVERGENCE_ORDER; ++n) {
      if (n > 0) {
        for (int d = 0; d < 3; ++d) {
          m_krnlPrototype.kDivMTSub(d,n) = seissol::init::kDivMTSub::Values[seissol::tensor::kDivMTSub::index(d,n)];
        }
      }
      m_krnlPrototype.selectModes(n) = seissol::init::selectModes::Values[seissol::tensor::selectModes::index(n)];
    }
    for (int k = 0; k < NUMBER_OF_QUANTITIES; k++) {
      m_krnlPrototype.selectQuantity(k) = seissol::init::selectQuantity::Values[seissol::tensor::selectQuantity::index(k)];
      m_krnlPrototype.selectQuantity_G(k) = seissol::init::selectQuantity_G::Values[seissol::tensor::selectQuantity_G::index(k)];
      //m_krnlPrototype.selectQuantity_Z(k) = seissol::init::selectQuantity_Z::Values[seissol::tensor::selectQuantity_Z::index(k)];
    }
    m_krnlPrototype.timeInt = seissol::init::timeInt::Values;
    m_krnlPrototype.wHat = seissol::init::wHat::Values;
  }

  void prepareQ(real* QData) {
    //scale quantities to make it more realistic
    std::array<real, 13> factor = {{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto Q = init::Q::view::create(QData);
    std::srand(1234);
    for (int  q = 0; q < NUMBER_OF_QUANTITIES; q++) {
      for (int bf = 0; bf < NUMBER_OF_BASIS_FUNCTIONS; bf++) {
        Q(bf, q) = (real)std::rand()/RAND_MAX*factor.at(q);
      }
    }
  }

  void solveWithKernel(real stp[]) {
    real o_timeIntegrated[seissol::tensor::I::size()];
    real stpRhs[seissol::tensor::stpRhs::size()] __attribute__((aligned(PAGESIZE_STACK)));
    std::fill(std::begin(stpRhs), std::end(stpRhs), 0);

    seissol::kernel::stp krnl;
    prepareKernel(krnl);

    real starMatrices0[tensor::star::size(0)];
    real starMatrices1[tensor::star::size(1)];
    real starMatrices2[tensor::star::size(2)];
    real sourceMatrix[tensor::ET::size()];
    real zMatrix[NUMBER_OF_QUANTITIES][tensor::Zinv::size(0)];
    prepareModel(starMatrices0, starMatrices1, starMatrices2, sourceMatrix, zMatrix);

    krnl.star(0) = starMatrices0;
    krnl.star(1) = starMatrices1;
    krnl.star(2) = starMatrices2;
    krnl.G = sourceMatrix;
    for(size_t i = 0; i < NUMBER_OF_QUANTITIES; i++) {
      krnl.Zinv(i)  = zMatrix[i];
    }

    real QData[NUMBER_OF_QUANTITIES*NUMBER_OF_BASIS_FUNCTIONS];
    prepareQ(QData);

    krnl.Q = QData;
    krnl.I = o_timeIntegrated;
    krnl.timestep = dt;
    krnl.stp = stp;
    krnl.stpRhs = stpRhs;
    krnl.execute();
  }

  public:
    void testSTP() {
#if CONVERGENCE_ORDER < 3 || CONVERGENCE_ORDER > 6
      std::cout << "STP test for order <3 or >6 not available" << std::endl;
      TS_ASSERT(false);
      return;
#endif
      real stp[seissol::tensor::stp::size()] __attribute__((aligned(PAGESIZE_STACK))) = {};
      std::fill(std::begin(stp), std::end(stp), 0);
      solveWithKernel(stp);
      auto withKernel = init::stp::view::create(stp);

      unsigned index = 0;
      for (unsigned i = 0; i < NUMBER_OF_BASIS_FUNCTIONS; i++) {
        for (unsigned j = 0; j < NUMBER_OF_QUANTITIES; j++) {
          for (unsigned k = 0; k < CONVERGENCE_ORDER; k++) {
            auto relativeError = std::abs(withKernel(i,j,k) - expectedValues.at(CONVERGENCE_ORDER).at(index)) /
              std::abs(expectedValues.at(CONVERGENCE_ORDER).at(index));
            TS_ASSERT(relativeError < epsilon*100);
            index++;
          }
        }
      }
    }
};

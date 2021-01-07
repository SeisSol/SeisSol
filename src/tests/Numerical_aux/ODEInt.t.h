#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "Numerical_aux/ODEVector.h"
#include <Numerical_aux/ODEInt.h>

namespace seissol::unit_test {
  class ODEIntTestSuite;
}

class seissol::unit_test::ODEIntTestSuite : public CxxTest::TestSuite {
private:
  double eps = 10e-10;
public:
  ODEIntTestSuite() {};
 void testSimpleIntegration() {
   constexpr auto sizeSolution = 5;
   constexpr auto sizeIntegratedSolution = sizeSolution;

   alignas(ALIGNMENT) real curUSolutionIntegrated[sizeIntegratedSolution] = {};
   alignas(ALIGNMENT) real fEvalIntegrated[sizeIntegratedSolution] = {};
   alignas(ALIGNMENT) real stage1Integrated[sizeIntegratedSolution] = {};
   alignas(ALIGNMENT) real fEvalHeunIntegrated[sizeIntegratedSolution] = {};
   alignas(ALIGNMENT) real updateHeunIntegrated[sizeIntegratedSolution] = {};

   alignas(ALIGNMENT) real curUSolution[sizeSolution] = {};
   alignas(ALIGNMENT) real fEvalEta[sizeSolution] = {};
   alignas(ALIGNMENT) real stage1Eta[sizeSolution] = {};
   alignas(ALIGNMENT) real fEvalHeunEta[sizeSolution] = {};
   alignas(ALIGNMENT) real updateHeunEta[sizeSolution] = {};

   auto curU = ODEVector{{curUSolutionIntegrated, curUSolution}, {sizeIntegratedSolution, sizeSolution}};
   auto fEval = ODEVector{{fEvalIntegrated, fEvalEta}, {sizeIntegratedSolution, sizeSolution}};
   auto stage1 = ODEVector{{stage1Integrated, stage1Eta}, {sizeIntegratedSolution, sizeSolution}};
   auto fEvalHeun = ODEVector{{fEvalHeunIntegrated, fEvalHeunEta}, {sizeIntegratedSolution, sizeSolution}};
   auto updateHeun = ODEVector{{updateHeunIntegrated, updateHeunEta}, {sizeIntegratedSolution, sizeSolution}};

   // Setup ODE solver
   const auto timeSpan = ode::TimeSpan{0, 2};
   const double dt = (timeSpan.end - timeSpan.begin) / 1;
   auto odeSolverConfig = ode::ODESolverConfig(dt);
   odeSolverConfig.acceptableError = 1e-10;
   odeSolverConfig.minimumDt = 1e-12;

   auto solver = ode::ODESolver(fEval, stage1, fEvalHeun, updateHeun, odeSolverConfig);
   auto f = [&](ODEVector& du,
               ODEVector& u,
               double time) {
     // Unpack du
     auto [dUIntegratedStorage, dUIntegratedSize] = du.getSubvector(0);
     assert(dUIntegratedSize == sizeIntegratedSolution);
     auto [dUStorage, dUSize] = du.getSubvector(1);
     assert(dUSize == sizeSolution);

     // Unpack u
     auto [uIntegrated, uIntegratedSize] = u.getSubvector(0);
     assert(uIntegratedSize == sizeIntegratedSolution);
     auto [uStorage, uSize] = u.getSubvector(1);
     assert(uSize == sizeSolution);

     // f(x)' = 2f(x) and f(0) = 1
     // Solution: exp(2x)
     // Solution for int u dx is exp(2x)/2 - C (C=0.5)
     du[0] = u[1];
     du[1] = 2*u[1];
     for (int i = 0; i < sizeSolution; ++i) {
       dUIntegratedStorage[i] = uStorage[i];
       dUStorage[i] = 2 * uStorage[i];
     }

     //du.print();
     //std::cout << "\n\n";
   };
   for (int i = 0; i < sizeSolution; ++i) {
     curUSolution[i] = 1.0;
     curUSolutionIntegrated[i] = 0.0;
   }

   solver.solve(f, curU, timeSpan);
   const double uShould = std::exp(2 * timeSpan.end);
   const double uIntegratedShould = std::exp(2 * timeSpan.end)/2 - 0.5;
   for (int i = 0; i < sizeSolution; ++i) {
     TS_ASSERT_DELTA(curUSolution[i], uShould, eps);
     TS_ASSERT_DELTA(curUSolutionIntegrated[i], uIntegratedShould, eps);
   }

 }
 };


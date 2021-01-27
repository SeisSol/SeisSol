#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "Numerical_aux/ODEVector.h"
#include <Numerical_aux/ODEInt.h>

namespace seissol::unit_test {
  class ODEIntTestSuite;
}

class seissol::unit_test::ODEIntTestSuite : public CxxTest::TestSuite {
private:
  double eps = 10e-11;
public:
  ODEIntTestSuite() {};
 void testSimpleIntegration() {
   constexpr auto sizeSolution = 5;
   constexpr auto sizeIntegratedSolution = sizeSolution;

   alignas(ALIGNMENT) real curUSolutionIntegrated[sizeIntegratedSolution] = {};
   alignas(ALIGNMENT) real curUSolution[sizeSolution] = {};

   auto curU = ODEVector{{curUSolutionIntegrated, curUSolution}, {sizeIntegratedSolution, sizeSolution}};

   // Setup ODE solver
   const auto timeSpan = ode::TimeSpan{0, 2};
   const double dt = 0.0005;
   auto odeSolverConfig = ode::ODESolverConfig(dt);
   odeSolverConfig.acceptableError = 1e-10;
   odeSolverConfig.minimumDt = 1e-12;

   auto solver = ode::RungeKuttaODESolver({sizeIntegratedSolution, sizeSolution}, odeSolverConfig);
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
     std::cout << std::abs(curUSolutionIntegrated[i] - uIntegratedShould) << std::endl;
   }

 }
 };


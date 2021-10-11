#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "Numerical_aux/ODEVector.h"
#include <Numerical_aux/ODEInt.h>

namespace seissol::unit_test {
class ODEIntTestSuite;
}

class seissol::unit_test::ODEIntTestSuite : public CxxTest::TestSuite {
private:
#ifdef SINGLE_PRECISION
  const real eps = 10e-4;
#else
  const real eps = 10e-11;
#endif
public:
  ODEIntTestSuite() {};

  void testSimpleIntegration() {
    constexpr auto sizeSolution = 5;
    constexpr auto sizeIntegratedSolution = sizeSolution;

    alignas(ALIGNMENT) real curUSolutionIntegrated[sizeIntegratedSolution] = {};
    alignas(ALIGNMENT) real curUSolution[sizeSolution] = {};

    auto curU = ode::ODEVector{{curUSolutionIntegrated, curUSolution},
                               {sizeIntegratedSolution, sizeSolution}};

    // Setup ODE solver
    const auto timeSpan = ode::TimeSpan{0, 2};
    const double dt = 0.0007;
    auto odeSolverConfig = ode::ODESolverConfig(dt);

    auto solver = ode::RungeKuttaODESolver({sizeIntegratedSolution, sizeSolution}, odeSolverConfig);
    auto f = [&](ode::ODEVector& du,
                 ode::ODEVector& u,
                 double time) {
      // Unpack du
      auto[dUIntegratedStorage, dUIntegratedSize] = du.getSubvector(0);
      assert(dUIntegratedSize == sizeIntegratedSolution);
      auto[dUStorage, dUSize] = du.getSubvector(1);
      assert(dUSize == sizeSolution);

      // Unpack u
      auto[uIntegrated, uIntegratedSize] = u.getSubvector(0);
      assert(uIntegratedSize == sizeIntegratedSolution);
      auto[uStorage, uSize] = u.getSubvector(1);
      assert(uSize == sizeSolution);

      // f(x)' = 2f(x) and f(0) = 1
      // Solution: exp(2x)
      // Solution for int u dx is exp(2x)/2 - C (C=0.5)
      //du[0] = u[1];
      //du[1] = 2 * u[1];
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
    const double uIntegratedShould = std::exp(2 * timeSpan.end) / 2 - 0.5;
    for (int i = 0; i < sizeSolution; ++i) {
      TS_ASSERT_DELTA(curUSolution[i], uShould, eps);
      TS_ASSERT_DELTA(curUSolutionIntegrated[i], uIntegratedShould, eps);
    }
  }


  void testLotkaVoltera() {
    constexpr auto sizeSolution = 2;

    alignas(ALIGNMENT) real curUSolution[sizeSolution] = {};

    auto curU = ode::ODEVector{{curUSolution},
                               {sizeSolution}};

    // Setup ODE solver
    const auto timeSpan = ode::TimeSpan{0, 1};
    const double dt = 0.01;
    auto odeSolverConfig = ode::ODESolverConfig(dt);
    odeSolverConfig.solver = ode::RungeKuttaVariant::RK7_VernerMostEfficient;

    auto solver = ode::RungeKuttaODESolver({sizeSolution}, odeSolverConfig);
    auto parameters = std::array<real, 4>{1.5, 1.0, 3.0, 1.0};
    auto f = [&](ode::ODEVector& du,
        ode::ODEVector& u,
        double time) {
      // A simple Lotka-Volterra model
      // See: https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations

      // dx/dt = \alpha x - \beta x y
      // dy/dt = \delta x y - \gamma y
      auto& [alpha, beta, delta, gamma] = parameters;
      du[0] = alpha * u[0] - beta * u[0] * u[1];
      du[1] = delta * u[0] * u[1] - gamma * u[1];
    };

    // Initial conditions
    curU[0] = 10.0;
    curU[1] = 1.0;

    solver.solve(f, curU, timeSpan);

    // Compare with reference solution computed with Vern 7 from DifferentialsEquations.jl
    /*
    function lotka_volterra!(du,u,p,t)
        x, y = u
        du[1] = p[1] * x - p[2] * x * y
        du[2] = p[3] * x * u[2] - p[4] * y
    end

    u0 = [10, 1.0]
    tspan = (0.0,1.0)
    p = [1.5, 1.0, 3.0, 1.0]
    prob = ODEProblem(lotka_volterra!,u0,tspan,p)
    sol = solve(prob, DifferentialEquations.Vern7(), dt=0.01, adaptive=false)
     */
    const auto uShould = std::array<double, 2>{2.070817357298899e-8, 15.074149470779822};
    TS_ASSERT_DELTA(curU[0], uShould[0], eps);
    TS_ASSERT_DELTA(curU[1], uShould[1], eps);
  }
};


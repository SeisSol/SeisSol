#ifndef SEISSOL_ODEINT_H
#define SEISSOL_ODEINT_H

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <tuple>

#include "ODEVector.h"

namespace seissol::ode {

struct TimeSpan {
  double begin;
  double end;
};

struct ODESolverConfig {
  double initialDt;
  double minimumDt = 1e-14;
  double acceptableError = 1e-14;

  ODESolverConfig() = delete;
  explicit ODESolverConfig(double initialDt) : initialDt(initialDt) {};
};

class ODESolver {
public:
  ODESolver(ODEVector& fEval,
            ODEVector& stage1,
            ODEVector& fEvalHeun,
            ODEVector& updateHeun,
            ODESolverConfig config
            ) :
            fEval(fEval),
            stage1(stage1),
            fEvalHeun(fEvalHeun),
            updateHeun(updateHeun),
            config(std::move(config)) {}
private:
  ODEVector& fEval;
  ODEVector& stage1;
  ODEVector& fEvalHeun;
  ODEVector& updateHeun;
  ODESolverConfig config;

public:
  template <typename Func>
  void solve(Func f,
             ODEVector& curValue,
             TimeSpan timeSpan) {
    assert(timeSpan.begin <= timeSpan.end);
    double curTime = timeSpan.begin;
    //curValue = initialValue;
    double dt = config.initialDt;
    //std::cout << "Begin ODE solve" << std::endl;
    while (curTime < timeSpan.end) {
      if (dt < config.minimumDt) {
        throw std::runtime_error("ODE solver: time step size smaller than minimal acceptable. Aborting!");
      }
      const double adjustedDt = std::min(dt, timeSpan.end - timeSpan.begin);
      //std::cout << "\ncurValue\n";
      //curValue.print();
      //std::cout << "\nfEval\n";
      f(fEval, curValue, curTime);
      //fEval.print();

      //const auto stage1 = curValue + curDt * fEval;
      stage1 = fEval; // = f(x_n, t_n)
      stage1 *= adjustedDt; // \hat{x_n+1} = adjustedDt * f(x_n, t_n)
      stage1 += curValue; // \hat{x_n+1} = curValue + adjustedDt * f(x_n, t_n)

      // TODO(Lukas) Can be optimized a lot...
      f(fEvalHeun, stage1, curTime + adjustedDt);
      //const auto updateHeun = curValue + 0.5 * dt * (
          //fEval + f(stage1, curTime)
          //);
      updateHeun = fEvalHeun;  // x_n+1 = f(\hat[x_n+1}, t_n+1)
      updateHeun += fEval;  // x_n+1 = f(\hat[x_n+1}, t_n+1) + f(x_n, t_n)
      updateHeun *= 0.5 * adjustedDt; // x_n+1 = 1/2 * dt * (f(\hat[x_n+1}, t_n+1) + f(x_n, t_n))
      updateHeun += curValue; // x_n+1 = x_n + 1/2 * dt * (f(\hat[x_n+1}, t_n+1) + f(x_n, t_n))

      const auto errorEstimate = updateHeun.normDifferenceTo(stage1);
      const auto updateSizeHeun = updateHeun.normDifferenceTo(curValue);
      const auto updateSizeEuler = stage1.normDifferenceTo(curValue);
      const auto fVal = fEval.norm();

      if (errorEstimate > config.acceptableError) {
        std::cout << "Error estimate: " << errorEstimate << " new dt = " << dt/2 << std::endl;
        dt /= 2;
      } else {
        //curValue = updateHeun;
        curValue = updateHeun; // TODO(Lukas) Use Heun update!
        curTime += adjustedDt;
      }
    }
  }
};

} // namespace seissol::ode

#endif //SEISSOL_ODEINT_H

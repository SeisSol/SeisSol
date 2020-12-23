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
             ODEVector curValue,
             TimeSpan timeSpan) {
    assert(timeSpan.begin <= timeSpan.end);
    double curTime = timeSpan.begin;
    //curValue = initialValue;
    double dt = config.initialDt;
    while (curTime < timeSpan.end) {
      if (dt < config.minimumDt) {
        throw std::runtime_error("ODE solver: time step size smaller than minimal acceptable. Aborting!");
      }
      const double adjustedDt = std::min(dt, timeSpan.end - timeSpan.begin);
      //const auto fEval = f(curValue, curTime);
      f(fEval, curValue, curTime);
      stage1 = fEval;

      //const auto stage1 = curValue + curDt * fEval;
      stage1 = fEval;
      stage1 *= adjustedDt;
      stage1 += curValue;

      // TODO(Lukas) Can be optimized a lot...
      f(fEvalHeun, stage1, curTime + adjustedDt);
      //const auto updateHeun = curValue + 0.5 * dt * (
          //fEval + f(stage1, curTime)
          //);
      updateHeun = fEvalHeun;
      updateHeun += fEval;
      updateHeun *= 0.5 * adjustedDt;

      // const auto errorEstimate = std::abs(stage1 - updateHeun);
      const auto errorEstimate = updateHeun.normDifferenceTo(stage1);

      if (errorEstimate > config.acceptableError) {
        dt /= 2;
      } else {
        curValue = updateHeun;
        curTime += adjustedDt;
      }
    }
  }
};

} // namespace seissol::ode

#endif //SEISSOL_ODEINT_H

#ifndef SEISSOL_ODEINT_H
#define SEISSOL_ODEINT_H

#include "Eigen/Dense"

#include "ODEVector.h"
#include "Kernels/precision.hpp"

namespace seissol::ode {

struct TimeSpan {
  double begin;
  double end;
};

enum class RungeKuttaVariant {
  RK4,
  RK4_3_8,
  RK4_Ralston,
  RK6_Butcher_1,
  RK6_Butcher_2,
  RK7_VernerMostEfficient
};

struct ODESolverConfig {
  RungeKuttaVariant solver = RungeKuttaVariant::RK7_VernerMostEfficient;
  double initialDt;

  ODESolverConfig() = delete;

  explicit ODESolverConfig(double initialDt) : initialDt(initialDt) {};
};

int getNumberOfStages(RungeKuttaVariant variant);

void initializeRungeKuttaScheme(RungeKuttaVariant variant,
                                int& numberOfStages,
                                Eigen::MatrixXd& a,
                                Eigen::VectorXd& b,
                                Eigen::VectorXd& c);

class RungeKuttaODESolver {
private:
  ODESolverConfig config;
  int numberOfStages{};

  // Coefficients
  Eigen::MatrixXd a;
  Eigen::VectorXd b;
  Eigen::VectorXd c;

  // Temporary storage
  std::vector<ODEVector> stages;
  std::vector<std::vector<real>> storages{};
  ODEVector buffer;

public:
  RungeKuttaODESolver(const std::vector<std::size_t>& storageSizes,
                      ODESolverConfig config);

  void setConfig(ODESolverConfig newConfig);


  /*!
   * @tparam Func is a callable type (e.g. functor/lambda)
   * @param f is a function with arguments
        f(ODEVector& du, ODEVector& u, double evaluationTime).
        which sets the right hand side of the ODE in du
        and takes u, du and evaluationTime as input.
   * @param curValue is the current solution of the ODE
   * @param timeSpan is the time span in which the ODE should be solved.
   */
  template<typename Func>
  void solve(Func f, ODEVector& curValue, TimeSpan timeSpan) {
    assert(timeSpan.begin <= timeSpan.end);
    double curTime = timeSpan.begin;
    double dt = config.initialDt;
    while (curTime < timeSpan.end) {
      const double adjustedDt = std::min(dt, timeSpan.end - curTime);

      for (auto i = 0U; i < stages.size(); ++i) {
        buffer = curValue;
        // j < i due to explict RK scheme
        for (auto j = 0U; j < i; ++j) {
          if (a(i,j) != 0.0) {
            const auto curWeight = a(i, j) * adjustedDt;
            buffer.weightedAddInplace(curWeight, stages[j]);
          }
        }

        const double tEval = curTime + c[i] * adjustedDt;
        f(stages[i], buffer, tEval);
      }

      for (auto i = 0; i < numberOfStages; ++i) {
        const auto curWeight = b[i] * adjustedDt;
        curValue.weightedAddInplace(curWeight, stages[i]);
      }
      curTime += adjustedDt;
    }
  }
};

} // namespace seissol::ode

#endif //SEISSOL_ODEINT_H

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
  RK6_Butcher_2
};

struct ODESolverConfig {
  RungeKuttaVariant solver = RungeKuttaVariant::RK6_Butcher_2;
  double initialDt;

  ODESolverConfig() = delete;

  explicit ODESolverConfig(double initialDt) : initialDt(initialDt) {};
};

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
          const auto curWeight = a(i, j) * adjustedDt;
          buffer.weightedAddInplace(curWeight, stages[j]);
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

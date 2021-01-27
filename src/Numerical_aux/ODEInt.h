#ifndef SEISSOL_ODEINT_H
#define SEISSOL_ODEINT_H

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <tuple>

#include "Eigen/Dense"

#include "ODEVector.h"

namespace seissol::ode {

struct TimeSpan {
  double begin;
  double end;
};

enum class RungeKuttaVariant {
  RK4,
  RK4_3_8,
  RK4_Ralston,
};

struct ODESolverConfig {
  RungeKuttaVariant solver = RungeKuttaVariant::RK4_Ralston;
  double initialDt;
  double minimumDt = 1e-14;

  // Only if adaptive (currently unsupported)
  double acceptableError = 1e-14;

  ODESolverConfig() = delete;
  explicit ODESolverConfig(double initialDt) : initialDt(initialDt) {};
};

inline void initializeRungeKuttaScheme(RungeKuttaVariant variant,
                                int& numberOfStages,
                                Eigen::MatrixXd& a,
                                Eigen::VectorXd& b,
                                Eigen::VectorXd& c) {
  std::unordered_map<RungeKuttaVariant, int> variantToNumberOfStages = {
      {RungeKuttaVariant::RK4, 4},
      {RungeKuttaVariant::RK4_Ralston, 4},
      {RungeKuttaVariant:: RK4_3_8, 4}
  };
  numberOfStages = variantToNumberOfStages[variant];

  // Initialize coefficients
  a = Eigen::MatrixXd(numberOfStages, numberOfStages);
  b = Eigen::VectorXd(numberOfStages);
  c = Eigen::VectorXd(numberOfStages);

  switch(variant) {
    case RungeKuttaVariant::RK4:
      // The classical RK4
      a(1,0) = 1.0/2.0;
      a(2,0) = 0.0;
      a(2,1) = 1.0/2.0;
      a(3,0) = 0.0;
      a(3,1) = 0.0;
      a(3,2) = 1.0;

      b(0) = 1.0/6.0;
      b(1) = 1.0/3.0;
      b(2) = 1.0/3.0;
      b(3) = 1.0/6.0;

      c(0) = 0.0;
      c(1) = 1.0/2.0;
      c(2) = 1.0/2.0;
      c(3) = 1.0;
      break;
    case RungeKuttaVariant::RK4_3_8:
      // The also classical 3/8 rule
      a(1,0) = 1.0/3.0;
      a(2,0) = -1.0/3.0;
      a(2,1) = 1.0;
      a(3,0) = 1.0;
      a(3,1) = -1.0;
      a(3,2) = 1.0;

      b(0) = 1.0/8.0;
      b(1) = 3.0/8.0;
      b(2) = 3.0/8.0;
      b(3) = 1.0/8.0;

      c(0) = 0.0;
      c(1) = 1.0/3.0;
      c(2) = 2.0/3.0;
      c(3) = 1.0;
      break;
    case RungeKuttaVariant::RK4_Ralston:
      // Ralston's RK4, minimized truncation error. Coeffs stolen from:
      // https://github.com/SciML/DiffEqDevTools.jl/blob/b5aca9330cd1a1b6ffbdbdf33a7ea037f7b53699/src/ode_tableaus.jl#L235

      a(1,0)= 4.0/10.0;
      a(2,0)= (-2889.0 + 1428.0 * std::sqrt(5.0)) / 1024.0;
      a(2,1)= (3785.0 - 1620.0 * std::sqrt(5.0)) / 1024.0;
      a(3,0)= (-3365.0 + 2094.0 * std::sqrt(5.0)) / 6040.0;
      a(3,1)= (-975.0 - 3046.0 * std::sqrt(5.0)) / 2552.0;
      a(3,2)= (467040.0 + 203968.0 * std::sqrt(5.0)) / 240845.0;

      b(0) = (263.0 + 24.0 * std::sqrt(5.0)) / 1812.0;
      b(1) = (125.0 - 1000.0 * std::sqrt(5.0)) / 3828.0;
      b(2) = 1024.0 * (3346.0 + 1623.0 * std::sqrt(5.0)) / 5924787.0;
      b(3) = (30.0 - 4.0 * std::sqrt(5.0)) / 123.0;

      c(0) = 0.0;
      c(1) = 4.0/10.0;
      c(2) = (14.0 - 3.0 * std::sqrt(5)) / 16.0;
      c(3) = 1.0;
      break;
  }
}

class RungeKuttaODESolver {
private:
  ODESolverConfig config;
  int numberOfStages;

  // Coefficients
  Eigen::MatrixXd a;
  Eigen::VectorXd b;
  Eigen::VectorXd c;

  // Temporary storage
  std::vector<ODEVector> stages;
  std::vector<std::vector<real>> storages;
  ODEVector buffer;

public:
  RungeKuttaODESolver(const std::vector<std::size_t>& storageSizes,
                      ODESolverConfig config)
                      : config(config) {
    initializeRungeKuttaScheme(config.solver, numberOfStages, a, b, c);
    // Initialize storages for stages
    stages.reserve(numberOfStages);
    storages.reserve((numberOfStages + 1) * storageSizes.size()); // +1 due to buffer
    auto curStoragePtrs = std::vector<real*>(storageSizes.size());
    for (auto i = 0; i < numberOfStages; ++i) {
      curStoragePtrs.clear();
      for (auto j = 0U; j < storageSizes.size(); ++j) {
        curStoragePtrs.push_back(storages.emplace_back(std::vector<real>(storageSizes[j])).data());
      }
      stages.push_back(ODEVector(curStoragePtrs, storageSizes));
    }

    // Initialize buffer
    curStoragePtrs.clear();
    for (auto j = 0U; j < storageSizes.size(); ++j) {
      curStoragePtrs.push_back(storages.emplace_back(std::vector<real>(storageSizes[j])).data());
    }
    buffer.updateStoragesAndSizes(curStoragePtrs, storageSizes);

  }

public:
  template <typename Func>
  void solve(Func f, ODEVector& curValue, TimeSpan timeSpan) {
    assert(timeSpan.begin <= timeSpan.end);
    double curTime = timeSpan.begin;
    double dt = config.initialDt;
    while (curTime < timeSpan.end) {
      if (dt < config.minimumDt) {
        throw std::runtime_error("ODE solver: time step size smaller than minimal acceptable. Aborting!");
      }
      const double adjustedDt = std::min(dt, timeSpan.end - curTime);

      for (auto i = 0U; i < stages.size(); ++i) {
        buffer = curValue;
        // j < i due to explict RK scheme
        for (auto j = 0U; j < i; ++j) {
          const auto curWeight = a(i,j) * adjustedDt;
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

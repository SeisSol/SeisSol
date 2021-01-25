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

struct ODESolverConfig {
  double initialDt;
  double minimumDt = 1e-14;
  double acceptableError = 1e-14;

  ODESolverConfig() = delete;
  explicit ODESolverConfig(double initialDt) : initialDt(initialDt) {};
};

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
                      : config(config), numberOfStages(4) {
    // Initialize coefficients
    a = Eigen::MatrixXd(numberOfStages, numberOfStages);
    b = Eigen::VectorXd(numberOfStages);
    c = Eigen::VectorXd(numberOfStages);

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

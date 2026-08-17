// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_NUMERICAL_ODEINT_H_
#define SEISSOL_SRC_NUMERICAL_ODEINT_H_

#include "Kernels/Precision.h"
#include "ODEVector.h"

#include <Eigen/Dense>
#include <cassert>

namespace seissol::ode {

struct TimeSpan {
  double begin;
  double end;
};

enum class RungeKuttaVariant {
  RK4,
  RK438,
  RK4Ralston,
  RK6Butcher1,
  RK6Butcher2,
  RK7VernerMostEfficient
};

struct ODESolverConfig {
  RungeKuttaVariant solver = RungeKuttaVariant::RK7VernerMostEfficient;
  double initialDt;

  ODESolverConfig() = delete;

  explicit ODESolverConfig(double initialDt) : initialDt(initialDt) {};
};

class RungeKuttaODESolver {
  private:
  ODESolverConfig config_;
  int numberOfStages_{};

  // Coefficients
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;

  // Temporary storage
  std::vector<ODEVector> stages_;
  std::vector<std::vector<real>> storages_;
  ODEVector buffer_;

  public:
  RungeKuttaODESolver(const std::vector<std::size_t>& storageSizes, ODESolverConfig config);

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
  template <typename Func>
  void solve(Func f, ODEVector& curValue, TimeSpan timeSpan) {
    assert(timeSpan.begin <= timeSpan.end);
    double curTime = timeSpan.begin;
    const double dt = config_.initialDt;
    while (curTime < timeSpan.end) {
      const double adjustedDt = std::min(dt, timeSpan.end - curTime);

      for (auto i = 0U; i < stages_.size(); ++i) {
        buffer_.copyFrom(curValue);
        // j < i due to explict RK scheme
        for (auto j = 0U; j < i; ++j) {
          if (a_[i * numberOfStages_ + j] != 0.0) {
            const auto curWeight = a_[i * numberOfStages_ + j] * adjustedDt;
            buffer_.weightedAddInplace(curWeight, stages_[j]);
          }
        }

        const double tEval = curTime + c_[i] * adjustedDt;
        f(stages_[i], buffer_, tEval);
      }

      for (auto i = 0; i < numberOfStages_; ++i) {
        const auto curWeight = b_[i] * adjustedDt;
        curValue.weightedAddInplace(curWeight, stages_[i]);
      }
      curTime += adjustedDt;
    }
  }
};

} // namespace seissol::ode

#endif // SEISSOL_SRC_NUMERICAL_ODEINT_H_

#include "ODEInt.h"
#include <algorithm>
#include <cassert>

namespace seissol::ode {
void initializeRungeKuttaScheme(RungeKuttaVariant variant,
                                int& numberOfStages,
                                Eigen::MatrixXd& a,
                                Eigen::VectorXd& b,
                                Eigen::VectorXd& c) {
  std::unordered_map<RungeKuttaVariant, int> variantToNumberOfStages = {
      {RungeKuttaVariant::RK4,           4},
      {RungeKuttaVariant::RK4_Ralston,   4},
      {RungeKuttaVariant::RK4_3_8,       4},
      {RungeKuttaVariant::RK6_Butcher_1, 7},
      {RungeKuttaVariant::RK6_Butcher_2, 7}
  };
  numberOfStages = variantToNumberOfStages[variant];

  // Initialize coefficients
  a = Eigen::MatrixXd(numberOfStages, numberOfStages);
  b = Eigen::VectorXd(numberOfStages);
  c = Eigen::VectorXd(numberOfStages);

  switch (variant) {
    case RungeKuttaVariant::RK4:
      // The classical RK4
      a(1, 0) = 1.0 / 2.0;
      a(2, 0) = 0.0;
      a(2, 1) = 1.0 / 2.0;
      a(3, 0) = 0.0;
      a(3, 1) = 0.0;
      a(3, 2) = 1.0;

      b(0) = 1.0 / 6.0;
      b(1) = 1.0 / 3.0;
      b(2) = 1.0 / 3.0;
      b(3) = 1.0 / 6.0;

      c(0) = 0.0;
      c(1) = 1.0 / 2.0;
      c(2) = 1.0 / 2.0;
      c(3) = 1.0;
      break;
    case RungeKuttaVariant::RK4_3_8:
      // The also classical 3/8 rule
      a(1, 0) = 1.0 / 3.0;
      a(2, 0) = -1.0 / 3.0;
      a(2, 1) = 1.0;
      a(3, 0) = 1.0;
      a(3, 1) = -1.0;
      a(3, 2) = 1.0;

      b(0) = 1.0 / 8.0;
      b(1) = 3.0 / 8.0;
      b(2) = 3.0 / 8.0;
      b(3) = 1.0 / 8.0;

      c(0) = 0.0;
      c(1) = 1.0 / 3.0;
      c(2) = 2.0 / 3.0;
      c(3) = 1.0;
      break;
    case RungeKuttaVariant::RK4_Ralston:
      // Ralston's RK4, minimized truncation error. Coeffs stolen from:
      // https://github.com/SciML/DiffEqDevTools.jl/blob/b5aca9330cd1a1b6ffbdbdf33a7ea037f7b53699/src/ode_tableaus.jl#L235

      a(1, 0) = 4.0 / 10.0;
      a(2, 0) = (-2889.0 + 1428.0 * std::sqrt(5.0)) / 1024.0;
      a(2, 1) = (3785.0 - 1620.0 * std::sqrt(5.0)) / 1024.0;
      a(3, 0) = (-3365.0 + 2094.0 * std::sqrt(5.0)) / 6040.0;
      a(3, 1) = (-975.0 - 3046.0 * std::sqrt(5.0)) / 2552.0;
      a(3, 2) = (467040.0 + 203968.0 * std::sqrt(5.0)) / 240845.0;

      b(0) = (263.0 + 24.0 * std::sqrt(5.0)) / 1812.0;
      b(1) = (125.0 - 1000.0 * std::sqrt(5.0)) / 3828.0;
      b(2) = 1024.0 * (3346.0 + 1623.0 * std::sqrt(5.0)) / 5924787.0;
      b(3) = (30.0 - 4.0 * std::sqrt(5.0)) / 123.0;

      c(0) = 0.0;
      c(1) = 4.0 / 10.0;
      c(2) = (14.0 - 3.0 * std::sqrt(5)) / 16.0;
      c(3) = 1.0;
      break;
    case RungeKuttaVariant::RK6_Butcher_1:
      // Coeffs from:
      // https://github.com/SciML/DiffEqDevTools.jl/blob/b5aca9330cd1a1b6ffbdbdf33a7ea037f7b53699/src/ode_tableaus.jl#L1266
      // (J.C. Butcher, 1964)
      a(1, 0) = 1.0 / 2.0 - 1 / 10.0 * std::sqrt(5.0);
      a(2, 0) = -1.0 / 10.0 * std::sqrt(5.0);
      a(2, 1) = 1.0 / 2.0 + 1.0 / 5.0 * std::sqrt(5.0);
      a(3, 0) = 7.0 / 20.0 * std::sqrt(5.0) - 3.0 / 4.0;
      a(3, 1) = 1.0 / 4.0 * std::sqrt(5.0) - 1.0 / 4.0;
      a(3, 2) = 3.0 / 2.0 - 7.0 / 10.0 * std::sqrt(5.0);
      a(4, 0) = 1.0 / 12.0 - 1.0 / 60.0 * std::sqrt(5.0);
      a(4, 1) = 0.0;
      a(4, 2) = 1.0 / 6.0;
      a(4, 3) = 7.0 / 60.0 * std::sqrt(5.0) + 1.0 / 4.0;
      a(5, 0) = 1.0 / 60.0 * std::sqrt(5.0) + 1.0 / 12.0;
      a(5, 1) = 0.0;
      a(5, 2) = 3.0 / 4.0 - 5.0 / 12.0 * std::sqrt(5.0);
      a(5, 3) = 1.0 / 6.0;
      a(5, 4) = -1.0 / 2.0 + 3.0 / 10.0 * std::sqrt(5.0);
      a(6, 0) = 1.0 / 6.0;
      a(6, 1) = 0;
      a(6, 2) = -55.0 / 12.0 + 25.0 / 12.0 * std::sqrt(5.0);
      a(6, 3) = -7.0 / 12.0 * std::sqrt(5.0) - 25.0 / 12.0;
      a(6, 4) = 5.0 - 2.0 * std::sqrt(5.0);
      a(6, 5) = 5.0 / 2.0 + 1.0 / 2.0 * std::sqrt(5.0);

      b(0) = 1.0 / 12.0;
      b(1) = 0.0;
      b(2) = 0.0;
      b(3) = 0.0;
      b(4) = 5.0 / 12.0;
      b(5) = 5.0 / 12.0;
      b(6) = 1.0 / 12.0;

      c(0) = 0.0;
      c(1) = 1.0 / 2.0 - 1.0 / 10.0 * std::sqrt(5.0);
      c(2) = 1.0 / 2.0 + 1.0 / 10.0 * std::sqrt(5.0);
      c(3) = 1.0 / 2.0 - 1.0 / 10.0 * std::sqrt(5.0);
      c(4) = 1.0 / 2.0 + 1.0 / 10.0 * std::sqrt(5.0);
      c(5) = 1.0 / 2.0 - 1.0 / 10.0 * std::sqrt(5.0);
      c(6) = 1.0;
      break;
    case RungeKuttaVariant::RK6_Butcher_2:
      // Coeffs from:
      // https://github.com/SciML/DiffEqDevTools.jl/blob/b5aca9330cd1a1b6ffbdbdf33a7ea037f7b53699/src/ode_tableaus.jl#L1320
      // (J.C. Butcher, 1964)
      a(1, 0) = 1.0 / 3.0;
      a(2, 0) = 0.0;
      a(2, 1) = 2.0 / 3.0;
      a(3, 0) = 1.0 / 12.0;
      a(3, 1) = 1.0 / 3.0;
      a(3, 2) = -1.0 / 12.0;
      a(4, 0) = -1.0 / 16.0;
      a(4, 1) = 9.0 / 8.0;
      a(4, 2) = -3.0 / 16.0;
      a(4, 3) = -3.0 / 8.0;
      a(5, 0) = 0.0;
      a(5, 1) = 9.0 / 8.0;
      a(5, 2) = -3.0 / 8.0;
      a(5, 3) = -3.0 / 4.0;
      a(5, 4) = 1.0 / 2.0;
      a(6, 0) = 9.0 / 44.0;
      a(6, 1) = -9.0 / 11.0;
      a(6, 2) = 63.0 / 44.0;
      a(6, 3) = 18.0 / 11.0;
      a(6, 4) = 0.0;
      a(6, 5) = -16.0 / 11.0;

      b(0) = 11.0 / 120.0;
      b(1) = 0.0;
      b(2) = 27.0 / 40.0;
      b(3) = 27.0 / 40.0;
      b(4) = -4.0 / 15.0;
      b(5) = -4.0 / 15.0;
      b(6) = 11.0 / 120.0;

      c(0) = 0.0;
      c(1) = 1.0 / 3.0;
      c(2) = 2.0 / 3.0;
      c(3) = 1.0 / 3.0;
      c(4) = 1.0 / 2.0;
      c(5) = 1.0 / 2.0;
      c(6) = 1.0;
      break;
  }
}

seissol::ode::RungeKuttaODESolver::RungeKuttaODESolver(const std::vector<std::size_t>& storageSizes,
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

} // namespace seissol::ode
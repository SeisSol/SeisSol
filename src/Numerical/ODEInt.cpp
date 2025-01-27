// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ODEInt.h"
#include <Kernels/Precision.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace seissol::ode {

int getNumberOfStages(RungeKuttaVariant variant) {
  std::unordered_map<RungeKuttaVariant, int> variantToNumberOfStages = {
      {RungeKuttaVariant::RK4, 4},
      {RungeKuttaVariant::RK4Ralston, 4},
      {RungeKuttaVariant::RK438, 4},
      {RungeKuttaVariant::RK6Butcher1, 7},
      {RungeKuttaVariant::RK6Butcher2, 7},
      {RungeKuttaVariant::RK7VernerMostEfficient, 9}};
  return variantToNumberOfStages[variant];
}
void initializeRungeKuttaScheme(RungeKuttaVariant variant,
                                int& numberOfStages,
                                Eigen::MatrixXd& a,
                                Eigen::VectorXd& b,
                                Eigen::VectorXd& c) {
  numberOfStages = getNumberOfStages(variant);

  // Initialize coefficients
  a = Eigen::MatrixXd(numberOfStages, numberOfStages);
  b = Eigen::VectorXd(numberOfStages);
  c = Eigen::VectorXd(numberOfStages);

  a.setZero();
  b.setZero();
  c.setZero();

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
  case RungeKuttaVariant::RK438:
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
  case RungeKuttaVariant::RK4Ralston:
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
  case RungeKuttaVariant::RK6Butcher1:
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
  case RungeKuttaVariant::RK6Butcher2:
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
  case RungeKuttaVariant::RK7VernerMostEfficient:
    // From Jim Verner:
    // http://people.math.sfu.ca/~jverner/
    // Here without order 6 estimate and without interpolation
    a(1, 0) = 0.005;
    a(2, 0) = -1.07679012345679;
    a(3, 0) = 0.04083333333333333;
    a(4, 0) = 0.6389139236255726;
    a(5, 0) = -2.6615773750187572;
    a(6, 0) = 6.067741434696772;
    a(7, 0) = 12.054670076253203;
    a(8, 0) = 10.138146522881808;
    a(2, 1) = 1.185679012345679;
    a(3, 2) = 0.1225;
    a(4, 2) = -2.455672638223657;
    a(5, 2) = 10.804513886456137;
    a(6, 2) = -24.711273635911088;
    a(7, 2) = -49.75478495046899;
    a(8, 2) = -42.6411360317175;
    a(4, 3) = 2.272258714598084;
    a(5, 3) = -8.3539146573962;
    a(6, 3) = 20.427517930788895;
    a(7, 3) = 41.142888638604674;
    a(8, 3) = 35.76384003992257;
    a(5, 4) = 0.820487594956657;
    a(6, 4) = -1.9061579788166472;
    a(7, 4) = -4.461760149974004;
    a(8, 4) = -4.3480228403929075;
    a(6, 5) = 1.006172249242068;
    a(7, 5) = 2.042334822239175;
    a(8, 5) = 2.0098622683770357;
    a(7, 6) = -0.09834843665406107;
    a(8, 6) = 0.3487490460338272;
    a(8, 7) = -0.27143900510483127;

    b(0) = 0.04715561848627222;
    b(3) = 0.25750564298434153;
    b(4) = 0.2621665397741262;
    b(5) = 0.15216092656738558;
    b(6) = 0.49399691700324844;
    b(7) = -0.29430311714032503;
    b(8) = 0.0813174723249511;

    c(1) = 0.005;
    c(2) = 0.10888888888888888;
    c(3) = 0.16333333333333333;
    c(4) = 0.4555;
    c(5) = 0.6095094489978381;
    c(6) = 0.884;
    c(7) = 0.925;
    c(8) = 1.0;
    break;
  }
}

RungeKuttaODESolver::RungeKuttaODESolver(const std::vector<std::size_t>& storageSizes,
                                         ODESolverConfig config)
    : config(config) {
  initializeRungeKuttaScheme(config.solver, numberOfStages, a, b, c);
  // Initialize storages for stages
  stages.reserve(numberOfStages);
  storages.reserve((numberOfStages + 1) * storageSizes.size()); // +1 due to buffer
  auto curStoragePtrs = std::vector<real*>(storageSizes.size());
  for (auto i = 0; i < numberOfStages; ++i) {
    curStoragePtrs.clear();
    for (const unsigned long storageSize : storageSizes) {
      curStoragePtrs.push_back(storages.emplace_back(storageSize).data());
    }
    stages.emplace_back(curStoragePtrs, storageSizes);
  }

  // Initialize buffer
  curStoragePtrs.clear();
  for (const unsigned long storageSize : storageSizes) {
    curStoragePtrs.push_back(storages.emplace_back(storageSize).data());
  }
  buffer.updateStoragesAndSizes(curStoragePtrs, storageSizes);
}

void RungeKuttaODESolver::setConfig(ODESolverConfig newConfig) { config = newConfig; }

} // namespace seissol::ode

// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Functions.h"

#include <array>
#include <cmath>
#include <cstdint>

namespace seissol::functions {

uint64_t rangeProduct(uint64_t from, uint64_t to) {
  uint64_t product = 1;
  for (; from <= to; ++from) {
    product *= from;
  }
  return product;
}

double JacobiP(unsigned n, unsigned a, unsigned b, double x) {
  if (n == 0) {
    return 1.0;
  }
  double pm2 = 0.0;
  double pm1 = 1.0;
  double pm = 0.5 * a - 0.5 * b + (1.0 + 0.5 * (a + b)) * x;
  const double a2B2 = static_cast<double>(a * a) - static_cast<double>(b * b);
  for (unsigned m = 2; m <= n; ++m) {
    pm2 = pm1;
    pm1 = pm;
    pm = ((2.0 * m + a + b - 1.0) * (a2B2 + (2.0 * m + a + b) * (2.0 * m + a + b - 2.0) * x) * pm1 -
          2.0 * (m + a - 1.0) * (m + b - 1.0) * (2.0 * m + a + b) * pm2) /
         (2.0 * m * (m + a + b) * (2.0 * m + a + b - 2.0));
  }
  return pm;
}

double JacobiPDerivative(unsigned n, unsigned a, unsigned b, double x) {
  return (n == 0) ? 0.0 : 0.5 * (n + a + b + 1.0) * JacobiP(n - 1, a + 1, b + 1, x);
}

std::array<double, 5> SingularityFreeJacobiPFactors(unsigned m, unsigned a, unsigned b) {
  const double c0 = 2.0 * m + a + b;
  const double c1 = c0 - 1.0;
  const double c2 = static_cast<double>(a * a) - static_cast<double>(b * b);
  const double c3 = c0 * (c0 - 2.0);
  const double c4 = 2.0 * (m + a - 1.0) * (m + b - 1.0) * c0;
  const double c5 = 2.0 * m * (m + a + b) * (c0 - 2.0);
  return {c1, c2, c3, c4, c5};
}

double SingularityFreeJacobiPRecursion(
    double x, double y, const std::array<double, 5>& cm, double pm1, double pm2) {
  return (cm[0] * (cm[1] * y + cm[2] * x) * pm1 - cm[3] * y * y * pm2) / cm[4];
}

double SingularityFreeJacobiP(unsigned n, unsigned a, unsigned b, double x, double y) {
  if (n == 0) {
    return 1.0;
  }
  double pm2 = 0.0;
  double pm1 = 1.0;
  double pm = (0.5 * a - 0.5 * b) * y + (1.0 + 0.5 * (a + b)) * x;
  for (unsigned m = 2; m <= n; ++m) {
    pm2 = pm1;
    pm1 = pm;
    auto c = SingularityFreeJacobiPFactors(m, a, b);
    pm = SingularityFreeJacobiPRecursion(x, y, c, pm1, pm2);
  }
  return pm;
}

std::array<double, 3>
    SingularityFreeJacobiPAndDerivatives(unsigned n, unsigned a, unsigned b, double x, double y) {
  if (n == 0) {
    return {1.0, 0.0, 0.0};
  }
  double pm2 = 0.0;
  double ddxPm2 = 0.0;
  double ddyPm2 = 0.0;
  double pm1 = 1.0;
  double ddxPm1 = 0.0;
  double ddyPm1 = 0.0;
  double pm = SingularityFreeJacobiP(1, a, b, x, y);
  double ddxPm = 1.0 + 0.5 * (a + b);
  double ddyPm = 0.5 * (static_cast<double>(a) - static_cast<double>(b));
  for (unsigned m = 2; m <= n; ++m) {
    pm2 = pm1;
    pm1 = pm;
    ddxPm2 = ddxPm1;
    ddxPm1 = ddxPm;
    ddyPm2 = ddyPm1;
    ddyPm1 = ddyPm;
    auto c = SingularityFreeJacobiPFactors(m, a, b);
    pm = SingularityFreeJacobiPRecursion(x, y, c, pm1, pm2);
    ddxPm = (c[0] * (c[2] * pm1 + (c[1] * y + c[2] * x) * ddxPm1) - c[3] * y * y * ddxPm2) / c[4];
    ddyPm = (c[0] * (c[1] * pm1 + (c[1] * y + c[2] * x) * ddyPm1) -
             c[3] * (2.0 * y * pm2 + y * y * ddyPm2)) /
            c[4];
  }
  return {pm, ddxPm, ddyPm};
}

double TriDubinerP(const std::array<unsigned, 2>& i, const std::array<double, 2>& xi) {
  const double rNum = 2.0 * xi[0] - 1.0 + xi[1];
  const double s = 2.0 * xi[1] - 1.0;
  const double theta = 1.0 - xi[1];

  const double ti = SingularityFreeJacobiP(i[0], 0, 0, rNum, theta);
  const double tij = SingularityFreeJacobiP(i[1], 2 * i[0] + 1, 0, s, 1.0);

  return ti * tij;
}

std::array<double, 2> gradTriDubinerP(const std::array<unsigned, 2>& i,
                                      const std::array<double, 2>& xi) {
  const double rNum = 2.0 * xi[0] - 1.0 + xi[1];
  const double s = 2.0 * xi[1] - 1.0;
  const double theta = 1.0 - xi[1];

  auto ti = SingularityFreeJacobiPAndDerivatives(i[0], 0, 0, rNum, theta);
  auto tij = SingularityFreeJacobiPAndDerivatives(i[1], 2 * i[0] + 1, 0, s, 1.0);

  auto ddalpha = [&](double drNum, double dtheta, double dt) {
    return (ti[1] * drNum + ti[2] * dtheta) * tij[0] + ti[0] * tij[1] * dt;
  };

  return {ddalpha(2.0, 0.0, 0.0), ddalpha(1.0, -1.0, 2.0)};
}

double TetraDubinerP(const std::array<unsigned, 3>& i, const std::array<double, 3>& xi) {
  const double rNum = 2.0 * xi[0] - 1.0 + xi[1] + xi[2];
  const double sNum = 2.0 * xi[1] - 1.0 + xi[2];
  const double t = 2.0 * xi[2] - 1.0;
  const double sigmatheta = 1.0 - xi[1] - xi[2];
  const double theta = 1.0 - xi[2];

  const double ti = SingularityFreeJacobiP(i[0], 0, 0, rNum, sigmatheta);
  const double tij = SingularityFreeJacobiP(i[1], 2 * i[0] + 1, 0, sNum, theta);
  const double tijk = SingularityFreeJacobiP(i[2], 2 * i[0] + 2 * i[1] + 2, 0, t, 1.0);

  return ti * tij * tijk;
}

std::array<double, 3> gradTetraDubinerP(const std::array<unsigned, 3>& i,
                                        const std::array<double, 3>& xi) {
  const double rNum = 2.0 * xi[0] - 1.0 + xi[1] + xi[2];
  const double sNum = 2.0 * xi[1] - 1.0 + xi[2];
  const double t = 2.0 * xi[2] - 1.0;
  const double sigmatheta = 1.0 - xi[1] - xi[2];
  const double theta = 1.0 - xi[2];

  auto ti = SingularityFreeJacobiPAndDerivatives(i[0], 0, 0, rNum, sigmatheta);
  auto tij = SingularityFreeJacobiPAndDerivatives(i[1], 2 * i[0] + 1, 0, sNum, theta);
  auto tijk = SingularityFreeJacobiPAndDerivatives(i[2], 2 * i[0] + 2 * i[1] + 2, 0, t, 1.0);

  auto ddalpha = [&](double drNum, double dsigmatheta, double dsNum, double dtheta, double dt) {
    return (ti[1] * drNum + ti[2] * dsigmatheta) * tij[0] * tijk[0] +
           ti[0] * (tij[1] * dsNum + tij[2] * dtheta) * tijk[0] + ti[0] * tij[0] * (tijk[1] * dt);
  };

  return {ddalpha(2.0, 0.0, 0.0, 0.0, 0.0),
          ddalpha(1.0, -1.0, 2.0, 0.0, 0.0),
          ddalpha(1.0, -1.0, 1.0, -1.0, 2.0)};
}

template <>
double DubinerP<1U>(const std::array<unsigned, 1U>& i, const std::array<double, 1U>& xi) {
  return JacobiP(i[0], 0, 0, 2.0 * xi[0] - 1.0);
}
template <>
double DubinerP<2U>(const std::array<unsigned, 2U>& i, const std::array<double, 2U>& xi) {
  return TriDubinerP(i, xi);
}
template <>
double DubinerP<3U>(const std::array<unsigned, 3U>& i, const std::array<double, 3U>& xi) {
  return TetraDubinerP(i, xi);
}

template <>
std::array<double, 1U> gradDubinerP<1U>(const std::array<unsigned, 1U>& i,
                                        const std::array<double, 1U>& xi) {
  return {JacobiPDerivative(i[0], 0, 0, 2.0 * xi[0] - 1.0)};
}
template <>
std::array<double, 2U> gradDubinerP<2U>(const std::array<unsigned, 2U>& i,
                                        const std::array<double, 2U>& xi) {
  return gradTriDubinerP(i, xi);
}
template <>
std::array<double, 3U> gradDubinerP<3U>(const std::array<unsigned, 3U>& i,
                                        const std::array<double, 3U>& xi) {
  return gradTetraDubinerP(i, xi);
}

} // namespace seissol::functions

// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_NUMERICAL_TIMEBASIS_H_
#define SEISSOL_SRC_NUMERICAL_TIMEBASIS_H_

#include "Functions.h"
#include <cstddef>
#include <vector>
namespace seissol::numerical {

template <typename RealT>
class TimeBasis {
  public:
  virtual std::vector<RealT> derivative(double position, double timestep) const = 0;
  virtual std::vector<RealT> point(double position, double timestep) const = 0;
  virtual std::vector<RealT> integrate(double start, double end, double timestep) const = 0;

  std::vector<RealT> collocate(const std::vector<double>& points, double timestep) const {
    std::vector<RealT> data;
    for (const auto& point : points) {
      const auto local = this->point(point, timestep);
      data.insert(data.end(), local.begin(), local.end());
    }
    return data;
  }
};

template <typename RealT>
class MonomialBasis : public TimeBasis<RealT> {
  public:
  explicit MonomialBasis(std::size_t order) : order(order) {}

  std::vector<RealT> derivative(double position, double timestep) const override {
    std::vector<RealT> coeffs(order);
    coeffs[0] = 0;
    if (coeffs.size() > 1) {
      coeffs[1] = 1;
      double coeffCache = 1;
      for (std::size_t i = 1; i < order; ++i) {
        coeffCache *= position / i;
        coeffs[i + 1] = coeffCache;
      }
    }
    return coeffs;
  }

  std::vector<RealT> point(double position, double timestep) const override {
    std::vector<RealT> coeffs(order);
    coeffs[0] = 1;
    double coeffCache = 1;
    for (std::size_t i = 1; i < order; ++i) {
      coeffCache *= position / i;
      coeffs[i] = coeffCache;
    }
    return coeffs;
  }

  std::vector<RealT> integrate(double start, double end, double timestep) const override {
    std::vector<RealT> coeffs(order);
    double coeffStart = start;
    double coeffEnd = end;
    for (std::size_t i = 0; i < order; ++i) {
      coeffs[i] = coeffEnd - coeffStart;
      coeffStart *= start / (i + 2);
      coeffEnd *= end / (i + 2);
    }
    return coeffs;
  }

  private:
  std::size_t order;
};

template <typename RealT>
class LegendreBasis : public TimeBasis<RealT> {
  public:
  explicit LegendreBasis(std::size_t order) : order(order) {}

  std::vector<RealT> derivative(double position, double timestep) const override {
    const double tau = position / timestep;
    std::vector<RealT> data(order);
    for (std::size_t i = 0; i < order; ++i) {
      data[i] = seissol::functions::shiftedLegendre(i, tau, 1);
    }
    return data;
  }

  std::vector<RealT> point(double position, double timestep) const override {
    const double tau = position / timestep;
    std::vector<RealT> data(order);
    for (std::size_t i = 0; i < order; ++i) {
      data[i] = seissol::functions::shiftedLegendre(i, tau, 0);
    }
    return data;
  }

  std::vector<RealT> integrate(double start, double end, double timestep) const override {
    const double tauS = start / timestep;
    const double tauE = end / timestep;
    std::vector<RealT> data(order);
    for (std::size_t i = 0; i < order; ++i) {
      data[i] = seissol::functions::shiftedLegendre(i, tauE, -1) -
                seissol::functions::shiftedLegendre(i, tauS, -1);
    }
    return data;
  }

  private:
  std::size_t order;
};

} // namespace seissol::numerical
#endif // SEISSOL_SRC_NUMERICAL_TIMEBASIS_H_

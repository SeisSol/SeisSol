// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_TPCOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_TPCOMMON_H_

#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace seissol::dr::friction_law::tp {

/**
 * Logarithmic gridpoints as defined in Noda&Lapusta (14). These are the \f$\hat{l}\f$ for
 * ThermalPressurization.
 */
template <typename RealT = real>
class GridPoints {
  public:
  explicit GridPoints(std::size_t count) : values_(count) {
    for (std::size_t i = 0; i < count; ++i) {
      values_[i] = misc::TpMaxWaveNumber * std::exp(-misc::TpLogDz * (count - i - 1));
    }
  }

  const RealT& operator[](std::size_t i) const { return values_[i]; };
  [[nodiscard]] const std::vector<RealT>& data() const { return values_; }
  [[nodiscard]] std::size_t size() const { return values_.size(); }

  private:
  std::vector<RealT> values_;
};

/**
 * Inverse Fourier coefficients on the logarithmic grid.
 */
template <typename RealT = real>
class InverseFourierCoefficients {
  public:
  explicit InverseFourierCoefficients(std::size_t count) : values_(count) {
    const GridPoints<double> localGridPoints(count);

    for (std::size_t i = 1; i + 1 < count; ++i) {
      values_[i] = std::sqrt(2 / M_PI) * localGridPoints[i] * misc::TpLogDz;
    }
    values_[0] = std::sqrt(2 / M_PI) * localGridPoints[0] * (1 + misc::TpLogDz);
    if (count > 1) {
      values_[count - 1] = std::sqrt(2 / M_PI) * localGridPoints[count - 1] * 0.5 * misc::TpLogDz;
    }
  }

  const RealT& operator[](std::size_t i) const { return values_[i]; };
  [[nodiscard]] const std::vector<RealT>& data() const { return values_; }
  [[nodiscard]] std::size_t size() const { return values_.size(); }

  private:
  std::vector<RealT> values_;
};

/**
 * Stores the heat generation (without tauV) \f$\exp\left(\hat{l}^2/2\right) / \sqrt{2 \pi}\f$.
 */
template <typename RealT = real>
class GaussianHeatSource {
  public:
  explicit GaussianHeatSource(std::size_t count) : values_(count) {
    const GridPoints<double> localGridPoints(count);
    const double factor = 1 / std::sqrt(2.0 * M_PI);

    for (std::size_t i = 0; i < count; ++i) {
      const double heatGeneration = std::exp(-0.5 * misc::power<2>(localGridPoints[i]));
      values_[i] = factor * heatGeneration;
    }
  }

  const RealT& operator[](std::size_t i) const { return values_[i]; };
  [[nodiscard]] const std::vector<RealT>& data() const { return values_; }
  [[nodiscard]] std::size_t size() const { return values_.size(); }

  private:
  std::vector<RealT> values_;
};

} // namespace seissol::dr::friction_law::tp

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_TPCOMMON_H_

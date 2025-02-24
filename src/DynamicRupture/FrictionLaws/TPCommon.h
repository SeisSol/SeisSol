// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_TPCOMMON_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_TPCOMMON_H_

#include <array>
#include <cstddef>

#include "DynamicRupture/Misc.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::tp {

/**
 * Logarithmic gridpoints as defined in Noda&Lapusta (14). These are the \f$\hat{l}\f$ for
 * ThermalPressurization.
 */
template <size_t N, typename RealT = real>
class GridPoints {
  public:
  GridPoints() {
    for (size_t i = 0; i < N; ++i) {
      values[i] =
          misc::TpMaxWaveNumber * std::exp(-misc::TpLogDz * (misc::NumTpGridPoints - i - 1));
    }
  }
  const RealT& operator[](size_t i) const { return values[i]; };
  const std::array<RealT, N>& data() const { return values; }

  private:
  std::array<RealT, N> values;
};

/**
 * Inverse Fourier coefficients on the logarithmic grid.
 */
template <size_t N, typename RealT = real>
class InverseFourierCoefficients {
  public:
  constexpr InverseFourierCoefficients() {
    const GridPoints<N, double> localGridPoints;

    for (size_t i = 1; i < N - 1; ++i) {
      values[i] = std::sqrt(2 / M_PI) * localGridPoints[i] * misc::TpLogDz;
    }
    values[0] = std::sqrt(2 / M_PI) * localGridPoints[0] * (1 + misc::TpLogDz);
    values[N - 1] = std::sqrt(2 / M_PI) * localGridPoints[N - 1] * 0.5 * misc::TpLogDz;
  }
  const RealT& operator[](size_t i) const { return values[i]; };
  const std::array<RealT, N>& data() const { return values; }

  private:
  std::array<RealT, N> values;
};

/**
 * Stores the heat generation (without tauV) \f$\exp\left(\hat{l}^2/2\right) / \sqrt{2 \pi}\f$.
 */
template <size_t N, typename RealT = real>
class GaussianHeatSource {
  public:
  constexpr GaussianHeatSource() {
    const GridPoints<N, double> localGridPoints;
    const double factor = 1 / std::sqrt(2.0 * M_PI);

    for (size_t i = 0; i < N; ++i) {
      const double heatGeneration = std::exp(-0.5 * misc::power<2>(localGridPoints[i]));
      values[i] = factor * heatGeneration;
    }
  }
  const RealT& operator[](size_t i) const { return values[i]; };
  const std::array<RealT, N>& data() const { return values; }

  private:
  std::array<RealT, N> values;
};

} // namespace seissol::dr::friction_law::tp

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_GPUIMPL_THERMALPRESSURIZATION_TPCOMMON_H_

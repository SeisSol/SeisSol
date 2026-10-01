// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ThermalPressurization.h"

#include "DynamicRupture/FrictionLaws/TPCommon.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {

void ThermalPressurization::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  temperature_ = layerData.var<LTSThermalPressurization::Temperature>();
  pressure_ = layerData.var<LTSThermalPressurization::Pressure>();
  theta_ = layerData.var<LTSThermalPressurization::Theta>();
  sigma_ = layerData.var<LTSThermalPressurization::Sigma>();
  halfWidthShearZone_ = layerData.var<LTSThermalPressurization::HalfWidthShearZone>();
  hydraulicDiffusivity_ = layerData.var<LTSThermalPressurization::HydraulicDiffusivity>();
}

namespace {
/// tau * V, the shear heating a point produces at the given state
inline real shearHeating(real mu, real normalStress, real slipRate) {
  return -mu * normalStress * slipRate;
}
} // namespace

void ThermalPressurization::prepareFluidPressure(real deltaT, std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
    real temperatureDiffusion = 0.0;
    real temperatureSlope = 0.0;
    real pressureDiffusion = 0.0;
    real pressureSlope = 0.0;

    const real lambdaPrime =
        drParameters_.undrainedTPResponse * drParameters_.thermalDiffusivity /
        (hydraulicDiffusivity_[ltsFace][pointIndex] - drParameters_.thermalDiffusivity);

    for (std::uint32_t tpGridPointIndex = 0; tpGridPointIndex < drParameters_.tpGridPoints;
         ++tpGridPointIndex) {
      // Gaussian shear zone in spectral domain, normalized by w
      // \hat{l} / w
      const real squaredNormalizedTpGrid =
          misc::power<2>(gridPoints_[tpGridPointIndex] / halfWidthShearZone_[ltsFace][pointIndex]);

      // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
      const real thetaTpGrid = drParameters_.thermalDiffusivity * squaredNormalizedTpGrid;
      const real sigmaTpGrid = hydraulicDiffusivity_[ltsFace][pointIndex] * squaredNormalizedTpGrid;
      const real preExpTheta = -thetaTpGrid * deltaT;
      const real preExpSigma = -sigmaTpGrid * deltaT;
      const real exp1mTheta = -std::expm1(preExpTheta);
      const real exp1mSigma = -std::expm1(preExpSigma);

      const std::size_t gridIndex = ltsFace * drParameters_.tpGridPoints + tpGridPointIndex;
      const real scaledInverseFourierCoefficient =
          inverseFourierCoefficients_[tpGridPointIndex] / halfWidthShearZone_[ltsFace][pointIndex];

      // The spectral state carries over as + F(t) exp(-A dt) in equation (10), which knows nothing
      // of this step's shear heating.
      temperatureDiffusion +=
          scaledInverseFourierCoefficient * theta_[gridIndex][pointIndex] * std::exp(preExpTheta);
      pressureDiffusion +=
          scaledInverseFourierCoefficient * sigma_[gridIndex][pointIndex] * std::exp(preExpSigma);

      // The generation B/A * (1 - exp(-A dt)) is linear in tau * V, so the whole dependence of the
      // step on the slip rate sits in one scalar and the grid does not have to be walked again.
      // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
      const real omegaSlope = heatSource_[tpGridPointIndex];
      temperatureSlope += scaledInverseFourierCoefficient * omegaSlope /
                          (drParameters_.heatCapacity * thetaTpGrid) * exp1mTheta;
      pressureSlope += scaledInverseFourierCoefficient * omegaSlope *
                       (drParameters_.undrainedTPResponse + lambdaPrime) /
                       (drParameters_.heatCapacity * sigmaTpGrid) * exp1mSigma;
    }

    // Update pore pressure change: sigma = pore pressure + lambda' * temperature
    temperatureOffset_[pointIndex] = temperatureDiffusion + drParameters_.initialTemperature;
    temperatureSlope_[pointIndex] = temperatureSlope;
    pressureOffset_[pointIndex] =
        -(pressureDiffusion - lambdaPrime * temperatureDiffusion) + drParameters_.initialPressure;
    pressureSlope_[pointIndex] = -(pressureSlope - lambdaPrime * temperatureSlope);
  }
}

void ThermalPressurization::applyShearHeating(
    const std::array<real, misc::NumPaddedPoints>& normalStress,
    const real (*mu)[misc::NumPaddedPoints],
    const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
    std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
    const real tauV = shearHeating(
        mu[ltsFace][pointIndex], normalStress[pointIndex], slipRateMagnitude[pointIndex]);
    temperature_[ltsFace][pointIndex] =
        temperatureOffset_[pointIndex] + temperatureSlope_[pointIndex] * tauV;
    pressure_[ltsFace][pointIndex] =
        pressureOffset_[pointIndex] + pressureSlope_[pointIndex] * tauV;
  }
}

void ThermalPressurization::finalizeFluidPressure(
    const std::array<real, misc::NumPaddedPoints>& normalStress,
    const real (*mu)[misc::NumPaddedPoints],
    const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
    real deltaT,
    std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
    const real tauV = shearHeating(
        mu[ltsFace][pointIndex], normalStress[pointIndex], slipRateMagnitude[pointIndex]);
    const real lambdaPrime =
        drParameters_.undrainedTPResponse * drParameters_.thermalDiffusivity /
        (hydraulicDiffusivity_[ltsFace][pointIndex] - drParameters_.thermalDiffusivity);

    for (std::uint32_t tpGridPointIndex = 0; tpGridPointIndex < drParameters_.tpGridPoints;
         ++tpGridPointIndex) {
      const real squaredNormalizedTpGrid =
          misc::power<2>(gridPoints_[tpGridPointIndex] / halfWidthShearZone_[ltsFace][pointIndex]);
      const real thetaTpGrid = drParameters_.thermalDiffusivity * squaredNormalizedTpGrid;
      const real sigmaTpGrid = hydraulicDiffusivity_[ltsFace][pointIndex] * squaredNormalizedTpGrid;
      const real preExpTheta = -thetaTpGrid * deltaT;
      const real preExpSigma = -sigmaTpGrid * deltaT;

      const real omega = tauV * heatSource_[tpGridPointIndex];
      const std::size_t gridIndex = ltsFace * drParameters_.tpGridPoints + tpGridPointIndex;
      theta_[gridIndex][pointIndex] =
          theta_[gridIndex][pointIndex] * std::exp(preExpTheta) +
          omega / (drParameters_.heatCapacity * thetaTpGrid) * -std::expm1(preExpTheta);
      sigma_[gridIndex][pointIndex] = sigma_[gridIndex][pointIndex] * std::exp(preExpSigma) +
                                      omega * (drParameters_.undrainedTPResponse + lambdaPrime) /
                                          (drParameters_.heatCapacity * sigmaTpGrid) *
                                          -std::expm1(preExpSigma);
    }

    // the coefficients of the preceding pass describe this very state, so the output follows from
    // them rather than from a second sum over the grid
    temperature_[ltsFace][pointIndex] =
        temperatureOffset_[pointIndex] + temperatureSlope_[pointIndex] * tauV;
    pressure_[ltsFace][pointIndex] =
        pressureOffset_[pointIndex] + pressureSlope_[pointIndex] * tauV;
  }
}

} // namespace seissol::dr::friction_law::cpu

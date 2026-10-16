// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ThermalPressurization.h"

#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Memory/Descriptor/DynamicRupture.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {

template <typename Cfg>
void ThermalPressurization<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  temperature_ = layerData.var<LTSThermalPressurization::Temperature>(Cfg());
  pressure_ = layerData.var<LTSThermalPressurization::Pressure>(Cfg());
  theta_ = layerData.var<LTSThermalPressurization::Theta>(Cfg());
  sigma_ = layerData.var<LTSThermalPressurization::Sigma>(Cfg());
  halfWidthShearZone_ = layerData.var<LTSThermalPressurization::HalfWidthShearZone>(Cfg());
  hydraulicDiffusivity_ = layerData.var<LTSThermalPressurization::HydraulicDiffusivity>(Cfg());
}

namespace {
/// tau * V, the shear heating a point produces at the given state
template <typename RealT>
RealT shearHeating(RealT mu, RealT normalStress, RealT slipRate) {
  return -mu * normalStress * slipRate;
}
} // namespace

template <typename Cfg>
void ThermalPressurization<Cfg>::prepareFluidPressure(real deltaT, std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
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

template <typename Cfg>
void ThermalPressurization<Cfg>::applyShearHeating(
    const std::array<real, misc::NumPaddedPoints<Cfg>>& normalStress,
    const real (*mu)[misc::NumPaddedPoints<Cfg>],
    const std::array<real, misc::NumPaddedPoints<Cfg>>& slipRateMagnitude,
    std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
    const real tauV = shearHeating(
        mu[ltsFace][pointIndex], normalStress[pointIndex], slipRateMagnitude[pointIndex]);
    temperature_[ltsFace][pointIndex] =
        temperatureOffset_[pointIndex] + temperatureSlope_[pointIndex] * tauV;
    pressure_[ltsFace][pointIndex] =
        pressureOffset_[pointIndex] + pressureSlope_[pointIndex] * tauV;
  }
}

template <typename Cfg>
void ThermalPressurization<Cfg>::finalizeFluidPressure(
    const std::array<real, misc::NumPaddedPoints<Cfg>>& normalStress,
    const real (*mu)[misc::NumPaddedPoints<Cfg>],
    const std::array<real, misc::NumPaddedPoints<Cfg>>& slipRateMagnitude,
    real deltaT,
    std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
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

      // Noda & Lapusta (2010) equation (10), F(t) exp(-A dt) + B/A (1 - exp(-A dt)), written as
      // F(t) + (B/A - F(t)) (1 - exp(-A dt)) with expm1 alone: a mode at its steady state B/A stays
      // there in any precision, rather than drifting by the mismatch of exp and expm1 in every
      // step.
      const real omega = tauV * heatSource_[tpGridPointIndex];
      const real thetaSteady = omega / (drParameters_.heatCapacity * thetaTpGrid);
      const real sigmaSteady = omega * (drParameters_.undrainedTPResponse + lambdaPrime) /
                               (drParameters_.heatCapacity * sigmaTpGrid);
      const std::size_t gridIndex = ltsFace * drParameters_.tpGridPoints + tpGridPointIndex;
      const real thetaOld = theta_[gridIndex][pointIndex];
      const real sigmaOld = sigma_[gridIndex][pointIndex];
      theta_[gridIndex][pointIndex] =
          thetaOld + (thetaSteady - thetaOld) * -std::expm1(preExpTheta);
      sigma_[gridIndex][pointIndex] =
          sigmaOld + (sigmaSteady - sigmaOld) * -std::expm1(preExpSigma);
    }

    // the coefficients of the preceding pass describe this very state, so the output follows from
    // them rather than from a second sum over the grid
    temperature_[ltsFace][pointIndex] =
        temperatureOffset_[pointIndex] + temperatureSlope_[pointIndex] * tauV;
    pressure_[ltsFace][pointIndex] =
        pressureOffset_[pointIndex] + pressureSlope_[pointIndex] * tauV;
  }
}

#define SEISSOL_INSTANTIATE(Cfg) template class ThermalPressurization<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law::cpu

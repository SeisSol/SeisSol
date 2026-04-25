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

static const tp::GridPoints<misc::NumTpGridPoints> TpGridPoints;
static const tp::InverseFourierCoefficients<misc::NumTpGridPoints> TpInverseFourierCoefficients;
static const tp::GaussianHeatSource<misc::NumTpGridPoints> HeatSource;

void ThermalPressurization::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  temperature_ = layerData.var<LTSThermalPressurization::Temperature>();
  pressure_ = layerData.var<LTSThermalPressurization::Pressure>();
  theta_ = layerData.var<LTSThermalPressurization::Theta>();
  sigma_ = layerData.var<LTSThermalPressurization::Sigma>();
  halfWidthShearZone_ = layerData.var<LTSThermalPressurization::HalfWidthShearZone>();
  hydraulicDiffusivity_ = layerData.var<LTSThermalPressurization::HydraulicDiffusivity>();
}

void ThermalPressurization::calcFluidPressure(
    const std::array<real, misc::NumPaddedPoints>& normalStress,
    const real (*mu)[misc::NumPaddedPoints],
    const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
    real deltaT,
    bool saveTPinLTS,
    std::size_t ltsFace) {
#pragma omp simd
  for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
    real temperatureUpdate = 0.0;
    real pressureUpdate = 0.0;

    const real faultStrength = -mu[ltsFace][pointIndex] * normalStress[pointIndex];
    const real tauV = faultStrength * slipRateMagnitude[pointIndex];
    const real lambdaPrime =
        drParameters_->undrainedTPResponse * drParameters_->thermalDiffusivity /
        (hydraulicDiffusivity_[ltsFace][pointIndex] - drParameters_->thermalDiffusivity);

    for (uint32_t tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
         ++tpGridPointIndex) {
      // Gaussian shear zone in spectral domain, normalized by w
      // \hat{l} / w
      const real squaredNormalizedTpGrid =
          misc::power<2>(TpGridPoints[tpGridPointIndex] / halfWidthShearZone_[ltsFace][pointIndex]);

      // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
      const real thetaTpGrid = drParameters_->thermalDiffusivity * squaredNormalizedTpGrid;
      const real sigmaTpGrid = hydraulicDiffusivity_[ltsFace][pointIndex] * squaredNormalizedTpGrid;
      const real preExpTheta = -thetaTpGrid * deltaT;
      const real preExpSigma = -sigmaTpGrid * deltaT;
      const real expTheta = std::exp(preExpTheta);
      const real expSigma = std::exp(preExpSigma);
      const real exp1mTheta = -std::expm1(preExpTheta);
      const real exp1mSigma = -std::expm1(preExpSigma);

      // Temperature and pressure diffusion in spectral domain over timestep
      // This is + F(t) exp(-A dt) in equation (10)
      const real thetaDiffusion = theta_[ltsFace][tpGridPointIndex][pointIndex] * expTheta;
      const real sigmaDiffusion = sigma_[ltsFace][tpGridPointIndex][pointIndex] * expSigma;

      // Heat generation during timestep
      // This is B/A * (1 - exp(-A dt)) in Noda & Lapusta (2010) equation (10)
      // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
      const real omega = tauV * HeatSource[tpGridPointIndex];
      const real thetaGeneration = omega / (drParameters_->heatCapacity * thetaTpGrid) * exp1mTheta;
      const real sigmaGeneration = omega * (drParameters_->undrainedTPResponse + lambdaPrime) /
                                   (drParameters_->heatCapacity * sigmaTpGrid) * exp1mSigma;

      // Sum both contributions up
      const auto thetaNew = thetaDiffusion + thetaGeneration;
      const auto sigmaNew = sigmaDiffusion + sigmaGeneration;

      // Recover temperature and altered pressure using inverse Fourier transformation from the new
      // contribution
      const real scaledInverseFourierCoefficient =
          TpInverseFourierCoefficients[tpGridPointIndex] / halfWidthShearZone_[ltsFace][pointIndex];
      temperatureUpdate += scaledInverseFourierCoefficient * thetaNew;
      pressureUpdate += scaledInverseFourierCoefficient * sigmaNew;

      if (saveTPinLTS) {
        theta_[ltsFace][tpGridPointIndex][pointIndex] = thetaNew;
        sigma_[ltsFace][tpGridPointIndex][pointIndex] = sigmaNew;
      }
    }
    // Update pore pressure change: sigma = pore pressure + lambda' * temperature
    pressureUpdate -= lambdaPrime * temperatureUpdate;

    // Temperature and pore pressure change at single GP on the fault + initial values
    temperature_[ltsFace][pointIndex] = temperatureUpdate + drParameters_->initialTemperature;
    pressure_[ltsFace][pointIndex] = -pressureUpdate + drParameters_->initialPressure;
  }
}

} // namespace seissol::dr::friction_law::cpu

// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ThermalPressurization.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include <DynamicRupture/FrictionLaws/TPCommon.h>
#include <algorithm>
#include <array>
#include <cmath>

namespace seissol::dr::friction_law::cpu {

static const tp::GridPoints<misc::NumTpGridPoints> TpGridPoints;
static const tp::InverseFourierCoefficients<misc::NumTpGridPoints> TpInverseFourierCoefficients;
static const tp::GaussianHeatSource<misc::NumTpGridPoints> HeatSource;

void ThermalPressurization::copyLtsTreeToLocal(
    seissol::initializer::Layer& layerData,
    const seissol::initializer::DynamicRupture* const dynRup,
    real fullUpdateTime) {
  const auto* concreteLts =
      dynamic_cast<const seissol::initializer::ThermalPressurization*>(dynRup);
  temperature = layerData.var(concreteLts->temperature);
  pressure = layerData.var(concreteLts->pressure);
  theta = layerData.var(concreteLts->theta);
  sigma = layerData.var(concreteLts->sigma);
  thetaTmpBuffer = layerData.var(concreteLts->thetaTmpBuffer);
  sigmaTmpBuffer = layerData.var(concreteLts->sigmaTmpBuffer);
  faultStrength = layerData.var(concreteLts->faultStrength);
  halfWidthShearZone = layerData.var(concreteLts->halfWidthShearZone);
  hydraulicDiffusivity = layerData.var(concreteLts->hydraulicDiffusivity);
}

void ThermalPressurization::calcFluidPressure(
    const std::array<real, misc::NumPaddedPoints>& normalStress,
    const real (*mu)[misc::NumPaddedPoints],
    const std::array<real, misc::NumPaddedPoints>& slipRateMagnitude,
    real deltaT,
    bool saveTPinLTS,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
    // compute fault strength
    faultStrength[ltsFace][pointIndex] = -mu[ltsFace][pointIndex] * normalStress[pointIndex];

    std::copy(&theta[ltsFace][pointIndex][0],
              &theta[ltsFace][pointIndex][misc::NumTpGridPoints],
              &thetaTmpBuffer[ltsFace][pointIndex][0]);
    std::copy(&sigma[ltsFace][pointIndex][0],
              &sigma[ltsFace][pointIndex][misc::NumTpGridPoints],
              &sigmaTmpBuffer[ltsFace][pointIndex][0]);

    // use Theta/Sigma from last timestep
    updateTemperatureAndPressure(
        slipRateMagnitude[pointIndex], deltaT, pointIndex, timeIndex, ltsFace);

    // copy back to LTS tree, if necessary
    if (saveTPinLTS) {
      std::copy(&thetaTmpBuffer[ltsFace][pointIndex][0],
                &thetaTmpBuffer[ltsFace][pointIndex][misc::NumTpGridPoints],
                &theta[ltsFace][pointIndex][0]);
      std::copy(&sigmaTmpBuffer[ltsFace][pointIndex][0],
                &sigmaTmpBuffer[ltsFace][pointIndex][misc::NumTpGridPoints],
                &sigma[ltsFace][pointIndex][0]);
    }
  }
}

void ThermalPressurization::updateTemperatureAndPressure(real slipRateMagnitude,
                                                         real deltaT,
                                                         unsigned int pointIndex,
                                                         unsigned int timeIndex,
                                                         unsigned int ltsFace) {
  real temperatureUpdate = 0.0;
  real pressureUpdate = 0.0;

  const real tauV = faultStrength[ltsFace][pointIndex] * slipRateMagnitude;
  const real lambdaPrime =
      drParameters->undrainedTPResponse * drParameters->thermalDiffusivity /
      (hydraulicDiffusivity[ltsFace][pointIndex] - drParameters->thermalDiffusivity);

#pragma omp simd
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::NumTpGridPoints;
       tpGridPointIndex++) {
    // Gaussian shear zone in spectral domain, normalized by w
    // \hat{l} / w
    const real squaredNormalizedTpGrid =
        misc::power<2>(TpGridPoints[tpGridPointIndex] / halfWidthShearZone[ltsFace][pointIndex]);

    // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
    const real thetaTpGrid = drParameters->thermalDiffusivity * squaredNormalizedTpGrid;
    const real sigmaTpGrid = hydraulicDiffusivity[ltsFace][pointIndex] * squaredNormalizedTpGrid;
    const real preExpTheta = -thetaTpGrid * deltaT;
    const real preExpSigma = -sigmaTpGrid * deltaT;
    const real expTheta = std::exp(preExpTheta);
    const real expSigma = std::exp(preExpSigma);
    const real exp1mTheta = -std::expm1(preExpTheta);
    const real exp1mSigma = -std::expm1(preExpSigma);

    // Temperature and pressure diffusion in spectral domain over timestep
    // This is + F(t) exp(-A dt) in equation (10)
    const real thetaDiffusion = thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expTheta;
    const real sigmaDiffusion = sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expSigma;

    // Heat generation during timestep
    // This is B/A * (1 - exp(-A dt)) in Noda & Lapusta (2010) equation (10)
    // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
    const real omega = tauV * HeatSource[tpGridPointIndex];
    const real thetaGeneration = omega / (drParameters->heatCapacity * thetaTpGrid) * exp1mTheta;
    const real sigmaGeneration = omega * (drParameters->undrainedTPResponse + lambdaPrime) /
                                 (drParameters->heatCapacity * sigmaTpGrid) * exp1mSigma;

    // Sum both contributions up
    thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = thetaDiffusion + thetaGeneration;
    sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = sigmaDiffusion + sigmaGeneration;

    // Recover temperature and altered pressure using inverse Fourier transformation from the new
    // contribution
    const real scaledInverseFourierCoefficient =
        TpInverseFourierCoefficients[tpGridPointIndex] / halfWidthShearZone[ltsFace][pointIndex];
    temperatureUpdate +=
        scaledInverseFourierCoefficient * thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex];
    pressureUpdate +=
        scaledInverseFourierCoefficient * sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex];
  }
  // Update pore pressure change: sigma = pore pressure + lambda' * temperature
  pressureUpdate -= lambdaPrime * temperatureUpdate;

  // Temperature and pore pressure change at single GP on the fault + initial values
  temperature[ltsFace][pointIndex] = temperatureUpdate + drParameters->initialTemperature;
  pressure[ltsFace][pointIndex] = -pressureUpdate + drParameters->initialPressure;
}

} // namespace seissol::dr::friction_law::cpu

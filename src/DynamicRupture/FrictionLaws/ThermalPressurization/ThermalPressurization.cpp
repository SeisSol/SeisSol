#include "ThermalPressurization.h"
#include <DynamicRupture/Misc.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/Layer.hpp>
#include <Kernels/precision.hpp>
#include <algorithm>
#include <array>
#include <cmath>

namespace seissol::dr::friction_law {

static const GridPoints<misc::numberOfTPGridPoints> tpGridPoints;
static const InverseFourierCoefficients<misc::numberOfTPGridPoints> tpInverseFourierCoefficients;
static const GaussianHeatSource<misc::numberOfTPGridPoints> heatSource;

void ThermalPressurization::copyLtsTreeToLocal(
    seissol::initializer::Layer& layerData,
    seissol::initializer::DynamicRupture const* const dynRup,
    real fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializer::LTSRateAndStateThermalPressurization const* const>(dynRup);
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
    std::array<real, misc::numPaddedPoints> const& normalStress,
    real const (*mu)[misc::numPaddedPoints],
    std::array<real, misc::numPaddedPoints> const& slipRateMagnitude,
    real deltaT,
    bool saveTPinLTS,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    // compute fault strength
    faultStrength[ltsFace][pointIndex] = -mu[ltsFace][pointIndex] * normalStress[pointIndex];

    std::copy(&theta[ltsFace][pointIndex][0],
              &theta[ltsFace][pointIndex][misc::numberOfTPGridPoints],
              &thetaTmpBuffer[ltsFace][pointIndex][0]);
    std::copy(&sigma[ltsFace][pointIndex][0],
              &sigma[ltsFace][pointIndex][misc::numberOfTPGridPoints],
              &sigmaTmpBuffer[ltsFace][pointIndex][0]);

    // use Theta/Sigma from last timestep
    updateTemperatureAndPressure(
        slipRateMagnitude[pointIndex], deltaT, pointIndex, timeIndex, ltsFace);

    // copy back to LTS tree, if necessary
    if (saveTPinLTS) {
      std::copy(&thetaTmpBuffer[ltsFace][pointIndex][0],
                &thetaTmpBuffer[ltsFace][pointIndex][misc::numberOfTPGridPoints],
                &theta[ltsFace][pointIndex][0]);
      std::copy(&sigmaTmpBuffer[ltsFace][pointIndex][0],
                &sigmaTmpBuffer[ltsFace][pointIndex][misc::numberOfTPGridPoints],
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
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
       tpGridPointIndex++) {
    // Gaussian shear zone in spectral domain, normalized by w
    // \hat{l} / w
    const real squaredNormalizedTPGrid =
        misc::power<2>(tpGridPoints[tpGridPointIndex] / halfWidthShearZone[ltsFace][pointIndex]);

    // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
    const real expTheta =
        std::exp(-drParameters->thermalDiffusivity * deltaT * squaredNormalizedTPGrid);
    const real expSigma =
        std::exp(-hydraulicDiffusivity[ltsFace][pointIndex] * deltaT * squaredNormalizedTPGrid);

    // Temperature and pressure diffusion in spectral domain over timestep
    // This is + F(t) exp(-A dt) in equation (10)
    const real thetaDiffusion = thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expTheta;
    const real sigmaDiffusion = sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expSigma;

    // Heat generation during timestep
    // This is B/A * (1 - exp(-A dt)) in Noda & Lapusta (2010) equation (10)
    // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
    const real omega = tauV * heatSource[tpGridPointIndex];
    const real thetaGeneration =
        omega /
        (drParameters->heatCapacity * squaredNormalizedTPGrid * drParameters->thermalDiffusivity) *
        (1.0 - expTheta);
    const real sigmaGeneration = omega * (drParameters->undrainedTPResponse + lambdaPrime) /
                                 (drParameters->heatCapacity * squaredNormalizedTPGrid *
                                  hydraulicDiffusivity[ltsFace][pointIndex]) *
                                 (1.0 - expSigma);

    // Sum both contributions up
    thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = thetaDiffusion + thetaGeneration;
    sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = sigmaDiffusion + sigmaGeneration;

    // Recover temperature and altered pressure using inverse Fourier transformation from the new
    // contribution
    const real scaledInverseFourierCoefficient =
        tpInverseFourierCoefficients[tpGridPointIndex] / halfWidthShearZone[ltsFace][pointIndex];
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

} // namespace seissol::dr::friction_law

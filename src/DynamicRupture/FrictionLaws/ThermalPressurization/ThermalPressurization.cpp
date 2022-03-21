#include "ThermalPressurization.h"

namespace seissol::dr::friction_law {

static const GridPoints<misc::numberOfTPGridPoints> tpGridPoints;
static const InverseFourierCoefficients<misc::numberOfTPGridPoints>
    tpInverseFourierCoefficients;

void ThermalPressurization::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                               seissol::initializers::DynamicRupture* dynRup,
                                               real fullUpdateTime) {
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurization*>(dynRup);
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

void ThermalPressurization::setInitialFluidPressure(unsigned int ltsFace) {
}

void ThermalPressurization::calcFluidPressure(
    std::array<real, misc::numPaddedPoints> const& normalStress,
    real (*mu)[misc::numPaddedPoints],
    std::array<real, misc::numPaddedPoints>& slipRateMagnitude,
    real deltaT,
    bool saveTmpInTP,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

    // compute fault strength
    faultStrength[ltsFace][pointIndex] =
        -mu[ltsFace][pointIndex] * normalStress[pointIndex];

    std::copy(&theta[ltsFace][pointIndex][0], &theta[ltsFace][pointIndex][misc::numberOfTPGridPoints], &thetaTmpBuffer[ltsFace][pointIndex][0]);
    std::copy(&sigma[ltsFace][pointIndex][0], &sigma[ltsFace][pointIndex][misc::numberOfTPGridPoints], &sigmaTmpBuffer[ltsFace][pointIndex][0]);

    // use Theta/Sigma from last call in this update, dt/2 and new SR from NS
    updateTemperatureAndPressure(
        slipRateMagnitude[pointIndex], deltaT, pointIndex, timeIndex, ltsFace);

    if (saveTmpInTP) {
      std::copy(&thetaTmpBuffer[ltsFace][pointIndex][0], &thetaTmpBuffer[ltsFace][pointIndex][misc::numberOfTPGridPoints], &theta[ltsFace][pointIndex][0]);
      std::copy(&sigmaTmpBuffer[ltsFace][pointIndex][0], &sigmaTmpBuffer[ltsFace][pointIndex][misc::numberOfTPGridPoints], &sigma[ltsFace][pointIndex][0]);
    }
  }
}

/**
 * We implement the time stepping algorithm, described in Noda&Lapusta (2010) eq. 10 for pressure and temperature in the frequency domain
 * @param slipRateMagnitude
 * @param deltaT
 * @param pointIndex
 * @param timeIndex
 * @param ltsFace
 */
void ThermalPressurization::updateTemperatureAndPressure(real slipRateMagnitude,
                                                         real deltaT,
                                                         unsigned int pointIndex,
                                                         unsigned int timeIndex,
                                                         unsigned int ltsFace) {
  real temperatureUpdate = 0.0;
  real pressureUpdate = 0.0;

  real tauV = faultStrength[ltsFace][pointIndex] * slipRateMagnitude;
  real lambdaPrime = drParameters.porePressureChange * drParameters.thermalDiffusivity /
                     (hydraulicDiffusivity[ltsFace][pointIndex] - drParameters.thermalDiffusivity);

#pragma omp simd
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
       tpGridPointIndex++) {
    // Gaussian shear zone in spectral domain, normalized by w
    real tmp = misc::power<2>(tpGridPoints.at(tpGridPointIndex) /
                              halfWidthShearZone[ltsFace][pointIndex]);

    // This is exp(-A dt) in equation (10)
    real expTheta = std::exp(-drParameters.thermalDiffusivity * deltaT * tmp);
    real expSigma = std::exp(-hydraulicDiffusivity[ltsFace][pointIndex] * deltaT * tmp);

    // Temperature and pressure diffusion in spectral domain over timestep
    // This is + F(t) exp(-A dt) in equation (10)
    real thetaDiffusion = thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expTheta;
    real sigmaDiffusion = sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expSigma;

    // Heat generation during timestep
    // This is B/A * (1 - exp(-A dt)) in equation (10)
    real omega = heatSource(tauV, tpGridPointIndex);
    real thetaGeneration = omega / (drParameters.heatCapacity * tmp * drParameters.thermalDiffusivity) * (1.0 - expTheta);
    real sigmaGeneration = omega * (drParameters.porePressureChange + lambdaPrime) / (drParameters.heatCapacity * tmp * hydraulicDiffusivity[ltsFace][pointIndex]) * (1.0 - expSigma);

    // Sum both contributions up
    thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = thetaDiffusion + thetaGeneration;
    sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = sigmaDiffusion + sigmaGeneration;

    // Recover temperature and pressure using inverse Fourier transformation from the new contribution
    real scaledInverseFourierCoefficient = tpInverseFourierCoefficients.at(tpGridPointIndex) / halfWidthShearZone[ltsFace][pointIndex];
    temperatureUpdate += scaledInverseFourierCoefficient * thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex];
    pressureUpdate += scaledInverseFourierCoefficient * sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex];
  }
  // Update pore pressure change: sigma = pore pressure + lambda' * temperature
  pressureUpdate = pressureUpdate - lambdaPrime * temperatureUpdate;

  // Temperature and pore pressure change at single GP on the fault + initial values
  temperature[ltsFace][pointIndex] = temperatureUpdate + drParameters.initialTemperature;
  pressure[ltsFace][pointIndex] = -pressureUpdate + drParameters.initialPressure;
}

/**
 * Implement Noda&Lapusta (2010) eq. (13)
 */
real ThermalPressurization::heatSource(double tauV, unsigned int tpGridPointIndex) {
  real factor = tauV / std::sqrt(2.0 * M_PI);
  real heatGeneration = std::exp(-0.5 * misc::power<2>(tpGridPoints.at(tpGridPointIndex)));

  return heatGeneration * factor;
}
} // namespace seissol::dr::friction_law

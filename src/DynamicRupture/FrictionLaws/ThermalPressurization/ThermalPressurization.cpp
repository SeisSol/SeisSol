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
    const FaultStresses& faultStresses,
    real (*initialStressInFaultCS)[misc::numPaddedPoints][6],
    real (*mu)[misc::numPaddedPoints],
    real (*slipRateMagnitude)[misc::numPaddedPoints],
    real deltaT,
    bool saveTmpInTP,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

    // compute fault strength
    auto normalStress = faultStresses.normalStress[timeIndex][pointIndex] +
                        initialStressInFaultCS[ltsFace][pointIndex][0] - pressure[ltsFace][pointIndex];
    faultStrength[ltsFace][pointIndex] =
        -mu[ltsFace][pointIndex] * std::min(static_cast<real>(0.0), normalStress);

    std::copy(&theta[ltsFace][pointIndex][0], &theta[ltsFace][pointIndex][misc::numberOfTPGridPoints], &thetaTmpBuffer[ltsFace][pointIndex][0]);
    std::copy(&sigma[ltsFace][pointIndex][0], &sigma[ltsFace][pointIndex][misc::numberOfTPGridPoints], &sigmaTmpBuffer[ltsFace][pointIndex][0]);

    //! use Theta/Sigma from last call in this update, dt/2 and new SR from NS
    updateTemperatureAndPressure(
        slipRateMagnitude[ltsFace][pointIndex], deltaT, pointIndex, timeIndex, ltsFace);

    if (saveTmpInTP) {
      std::copy(&thetaTmpBuffer[ltsFace][pointIndex][0], &thetaTmpBuffer[ltsFace][pointIndex][misc::numberOfTPGridPoints], &theta[ltsFace][pointIndex][0]);
      std::copy(&sigmaTmpBuffer[ltsFace][pointIndex][0], &sigmaTmpBuffer[ltsFace][pointIndex][misc::numberOfTPGridPoints], &sigma[ltsFace][pointIndex][0]);
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

  real tauV = faultStrength[ltsFace][pointIndex] * slipRateMagnitude;
  real lambdaPrime = drParameters.porePressureChange * drParameters.thermalDiffusivity /
                     (hydraulicDiffusivity[ltsFace][pointIndex] - drParameters.thermalDiffusivity);

#pragma omp simd
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
       tpGridPointIndex++) {
    // Gaussian shear zone in spectral domain, normalized by w
    real tmp = misc::power<2>(tpGridPoints.at(tpGridPointIndex) /
                              halfWidthShearZone[ltsFace][pointIndex]);

    // 1. Calculate diffusion of the field at previous timestep
    // temperature
    real thetaCurrent = thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * std::exp(-drParameters.thermalDiffusivity * deltaT * tmp);
    // pore pressuredeltaT + lambda'*temp
    real sigmaCurrent =
        sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * std::exp(-hydraulicDiffusivity[ltsFace][pointIndex] * deltaT * tmp);

    // 2. Add current contribution and get new temperature
    real omegaTheta = heatSource(tmp, drParameters.thermalDiffusivity, deltaT, tpGridPointIndex, timeIndex);
    thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = thetaCurrent + (tauV / drParameters.heatCapacity) * omegaTheta;

    real omegaSigma =
        heatSource(tmp, hydraulicDiffusivity[ltsFace][pointIndex], deltaT, tpGridPointIndex, timeIndex);
    sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = sigmaCurrent + ((drParameters.porePressureChange + lambdaPrime) * tauV) /
                                                                           (drParameters.heatCapacity) * omegaSigma;

    // 3. Recover temperature and pressure using inverse Fourier transformation with the calculated
    // fourier coefficients new contribution
    temperatureUpdate += (tpInverseFourierCoefficients.at(tpGridPointIndex) /
                          halfWidthShearZone[ltsFace][pointIndex]) *
                         thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex];
    pressureUpdate += (tpInverseFourierCoefficients.at(tpGridPointIndex) /
                       halfWidthShearZone[ltsFace][pointIndex]) *
                      sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex];
  }
  // Update pore pressure change (sigma = pore pressure + lambda'*temp)
  // In the BIEM code (Lapusta) they use T without initial value
  pressureUpdate = pressureUpdate - lambdaPrime * temperatureUpdate;

  // Temp and pore pressure change at single GP on the fault + initial values
  temperature[ltsFace][pointIndex] = temperatureUpdate + drParameters.initialTemperature;
  pressure[ltsFace][pointIndex] = -pressureUpdate + drParameters.initialPressure;
}

/**
 * Original function in spatial domain:
 * \f[\omega = \frac{1}{w \sqrt{2\pi}}\cdot
 * \exp\left(-\frac{1}{2}\cdot\left(\frac{z}{h}\right)^2\right) \f] Function in the wavenumber
 * domain including additional factors in front of the heat source function \f[ \omega =
 * \frac{1}{\alpha\cdot k^2 \cdot \sqrt{2\pi}}\cdot \exp\left(-\frac{1}{2}\cdot \left(k
 * h\right)^2\right)\cdot\left(1-\exp\left(-\alpha\cdot dt\cdot
 * \left(\frac{k}{h}\right)^2\right)\right) \f] inserting \f$\frac{k}{h}\f$ (scaled) for \f$k\f$
 * cancels out \f$h\f$fluid.
 * @param tmp \f$\left(k/h\right)^2\f$
 * @param alpha \f$\alpha\f$
 */
real ThermalPressurization::heatSource(
    real tmp, real alpha, real deltaT, unsigned int tpGridPointIndex, unsigned int timeIndex) {
  real value = 1.0 / (alpha * tmp * (sqrt(2.0 * M_PI))) *
               std::exp(-0.5 * misc::power<2>(tpGridPoints.at(tpGridPointIndex))) *
               (1.0 - exp(-alpha * deltaT * tmp));
  return value;
}
} // namespace seissol::dr::friction_law

#include "ThermalPressurization.h"

namespace seissol::dr::friction_law {

static constexpr GridPoints<misc::numberOfTPGridPoints> tpGridPoints;
static constexpr InverseFourierCoefficients<misc::numberOfTPGridPoints>
    tpInverseFourierCoefficients;

void ThermalPressurization::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                               seissol::initializers::DynamicRupture* dynRup,
                                               real fullUpdateTime) {
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurization*>(dynRup);
  temperature = layerData.var(concreteLts->temperature);
  pressure = layerData.var(concreteLts->pressure);
  tpTheta = layerData.var(concreteLts->tpTheta);
  tpSigma = layerData.var(concreteLts->tpSigma);
  tpHalfWidthShearZone = layerData.var(concreteLts->tpHalfWidthShearZone);
  alphaHy = layerData.var(concreteLts->alphaHy);
}

void ThermalPressurization::setInitialFluidPressure(unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    localPressure[pointIndex] = pressure[ltsFace][pointIndex];
  }
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
                        initialStressInFaultCS[ltsFace][pointIndex][0] - localPressure[pointIndex];
    faultStrength[pointIndex] =
        -mu[ltsFace][pointIndex] * std::min(static_cast<real>(0.0), normalStress);

    for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
         tpGridPointIndex++) {
      // save original values as they are overwritten updateTemperatureAndPressure
      thetaTmp[tpGridPointIndex] = tpTheta[ltsFace][pointIndex][tpGridPointIndex];
      sigmaTmp[tpGridPointIndex] = tpSigma[ltsFace][pointIndex][tpGridPointIndex];
    }
    //! use Theta/Sigma from last call in this update, dt/2 and new SR from NS
    updateTemperatureAndPressure(
        slipRateMagnitude[ltsFace][pointIndex], deltaT, pointIndex, timeIndex, ltsFace);

    localPressure[pointIndex] = pressure[ltsFace][pointIndex];
    if (saveTmpInTP) {
      for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
           tpGridPointIndex++) {
        tpTheta[ltsFace][pointIndex][tpGridPointIndex] = thetaTmp[tpGridPointIndex];
        tpSigma[ltsFace][pointIndex][tpGridPointIndex] = sigmaTmp[tpGridPointIndex];
      }
    }
  }
}

void ThermalPressurization::updateTemperatureAndPressure(real slipRateMagnitude,
                                                         real deltaT,
                                                         unsigned int pointIndex,
                                                         unsigned int timeIndex,
                                                         unsigned int ltsFace) {
  real localTemperature = 0.0;
  real localPressure = 0.0;

  real tauV = faultStrength[pointIndex] * slipRateMagnitude;
  real lambdaPrime = drParameters.tpLambda * drParameters.alphaTh /
                     (alphaHy[ltsFace][pointIndex] - drParameters.alphaTh);

#pragma omp simd
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
       tpGridPointIndex++) {
    // Gaussian shear zone in spectral domain, normalized by w
    real tmp = misc::power<2>(tpGridPoints.at(tpGridPointIndex) /
                              tpHalfWidthShearZone[ltsFace][pointIndex]);

    // 1. Calculate diffusion of the field at previous timestep
    // temperature
    real thetaCurrent = thetaTmp[tpGridPointIndex] * std::exp(-drParameters.alphaTh * deltaT * tmp);
    // pore pressure + lambda'*temp
    real sigmaCurrent =
        sigmaTmp[tpGridPointIndex] * std::exp(-alphaHy[ltsFace][pointIndex] * deltaT * tmp);

    // 2. Add current contribution and get new temperature
    real omegaTheta = heatSource(tmp, drParameters.alphaTh, deltaT, tpGridPointIndex, timeIndex);
    thetaTmp[tpGridPointIndex] = thetaCurrent + (tauV / drParameters.rhoC) * omegaTheta;

    real omegaSigma =
        heatSource(tmp, alphaHy[ltsFace][pointIndex], deltaT, tpGridPointIndex, timeIndex);
    sigmaTmp[tpGridPointIndex] = sigmaCurrent + ((drParameters.tpLambda + lambdaPrime) * tauV) /
                                                    (drParameters.rhoC) * omegaSigma;

    // 3. Recover temperature and pressure using inverse Fourier transformation with the calculated
    // fourier coefficients new contribution
    localTemperature += (tpInverseFourierCoefficients[tpGridPointIndex] /
                         tpHalfWidthShearZone[ltsFace][pointIndex]) *
                        thetaTmp[tpGridPointIndex];
    localPressure += (tpInverseFourierCoefficients[tpGridPointIndex] /
                      tpHalfWidthShearZone[ltsFace][pointIndex]) *
                     sigmaTmp[tpGridPointIndex];
  }
  // Update pore pressure change (sigma = pore pressure + lambda'*temp)
  // In the BIEM code (Lapusta) they use T without initial value
  localPressure = localPressure - lambdaPrime * localTemperature;

  // Temp and pore pressure change at single GP on the fault + initial values
  temperature[ltsFace][pointIndex] = localTemperature + drParameters.initialTemperature;
  pressure[ltsFace][pointIndex] = -localPressure + drParameters.initialPressure;
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
  return 1.0 / (alpha * tmp * (sqrt(2.0 * M_PI))) *
         std::exp(-0.5 * misc::power<2>(tpGridPoints.at(tpGridPointIndex))) *
         (1.0 - exp(-alpha * deltaT * tmp));
}
} // namespace seissol::dr::friction_law

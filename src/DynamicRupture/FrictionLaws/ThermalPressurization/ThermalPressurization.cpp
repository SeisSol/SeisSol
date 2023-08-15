#include "ThermalPressurization.h"

namespace seissol::dr::friction_law {

template <typename Config>
void ThermalPressurization<Config>::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    RealT fullUpdateTime) {
  auto* concreteLts = dynamic_cast<
      seissol::initializers::LTSRateAndStateThermalPressurization<Config> const* const>(dynRup);
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

template <typename Config>
void ThermalPressurization<Config>::calcFluidPressure(
    std::array<RealT, misc::numPaddedPoints<Config>> const& normalStress,
    RealT const (*mu)[misc::numPaddedPoints<Config>],
    std::array<RealT, misc::numPaddedPoints<Config>> const& slipRateMagnitude,
    RealT deltaT,
    bool saveTPinLTS,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints<Config>; pointIndex++) {
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

template <typename Config>
void ThermalPressurization<Config>::updateTemperatureAndPressure(RealT slipRateMagnitude,
                                                                 RealT deltaT,
                                                                 unsigned int pointIndex,
                                                                 unsigned int timeIndex,
                                                                 unsigned int ltsFace) {
  RealT temperatureUpdate = 0.0;
  RealT pressureUpdate = 0.0;

  const RealT tauV = faultStrength[ltsFace][pointIndex] * slipRateMagnitude;
  const RealT lambdaPrime =
      drParameters->undrainedTPResponse * drParameters->thermalDiffusivity /
      (hydraulicDiffusivity[ltsFace][pointIndex] - drParameters->thermalDiffusivity);

#pragma omp simd
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < misc::numberOfTPGridPoints;
       tpGridPointIndex++) {
    // Gaussian shear zone in spectral domain, normalized by w
    // \hat{l} / w
    const RealT squaredNormalizedTPGrid =
        misc::power<2>(tpGridPoints[tpGridPointIndex] / halfWidthShearZone[ltsFace][pointIndex]);

    // This is exp(-A dt) in Noda & Lapusta (2010) equation (10)
    const RealT expTheta =
        std::exp(-drParameters->thermalDiffusivity * deltaT * squaredNormalizedTPGrid);
    const RealT expSigma =
        std::exp(-hydraulicDiffusivity[ltsFace][pointIndex] * deltaT * squaredNormalizedTPGrid);

    // Temperature and pressure diffusion in spectral domain over timestep
    // This is + F(t) exp(-A dt) in equation (10)
    const RealT thetaDiffusion = thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expTheta;
    const RealT sigmaDiffusion = sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] * expSigma;

    // Heat generation during timestep
    // This is B/A * (1 - exp(-A dt)) in Noda & Lapusta (2010) equation (10)
    // heatSource stores \exp(-\hat{l}^2 / 2) / \sqrt{2 \pi}
    const RealT omega = tauV * heatSource[tpGridPointIndex];
    const RealT thetaGeneration =
        omega /
        (drParameters->heatCapacity * squaredNormalizedTPGrid * drParameters->thermalDiffusivity) *
        (1.0 - expTheta);
    const RealT sigmaGeneration = omega * (drParameters->undrainedTPResponse + lambdaPrime) /
                                  (drParameters->heatCapacity * squaredNormalizedTPGrid *
                                   hydraulicDiffusivity[ltsFace][pointIndex]) *
                                  (1.0 - expSigma);

    // Sum both contributions up
    thetaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = thetaDiffusion + thetaGeneration;
    sigmaTmpBuffer[ltsFace][pointIndex][tpGridPointIndex] = sigmaDiffusion + sigmaGeneration;

    // Recover temperature and altered pressure using inverse Fourier transformation from the new
    // contribution
    const RealT scaledInverseFourierCoefficient =
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

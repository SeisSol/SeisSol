#include "FastVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law {
void FastVelocityWeakeningLaw::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                                  seissol::initializers::DynamicRupture* dynRup,
                                                  real fullUpdateTime) {
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);

  averagedSlip = layerData.var(concreteLts->averagedSlip);
  srW = layerData.var(concreteLts->rsSrW);
}

real FastVelocityWeakeningLaw::updateStateVariable(unsigned int pointIndex,
                                                   unsigned int face,
                                                   real stateVarReference,
                                                   real timeIncrement,
                                                   real localSlipRate) {
  double muW = drParameters.muW;
  double localSrW = srW[face][pointIndex];
  double localA = a[face][pointIndex];
  double localSl0 = sl0[face][pointIndex];

  // low-velocity steady state friction coefficient
  real lowVelocityFriction =
      drParameters.rsF0 - (drParameters.rsB - localA) * log(localSlipRate / drParameters.rsSr0);
  real steadyStateFrictionCoefficient =
      muW + (lowVelocityFriction - muW) /
                std::pow(1.0 + misc::power<8>(localSlipRate / localSrW), 1.0 / 8.0);
  // For compiling reasons we write SINH(X)=(EXP(X)-EXP(-X))/2
  real steadyStateStateVariable = localA * log(2.0 * drParameters.rsSr0 / localSlipRate *
                                               (exp(steadyStateFrictionCoefficient / localA) -
                                                exp(-steadyStateFrictionCoefficient / localA)) /
                                               2.0);

  // exact integration of dSV/dt DGL, assuming constant V over integration step

  real exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
  real localStateVariable = steadyStateStateVariable * (1.0 - exp1) + exp1 * stateVarReference;
  assert(!(std::isnan(localStateVariable) && pointIndex < misc::numberOfBoundaryGaussPoints) &&
         "NaN detected");
  return localStateVariable;
}

real FastVelocityWeakeningLaw::updateMu(unsigned int ltsFace,
                                        unsigned int pointIndex,
                                        real localSlipRateMagnitude,
                                        real localStateVariable) {
  // mu = a * arcsinh ( V / (2*V_0) * exp (psi / a))
  real localA = a[ltsFace][pointIndex];
  // x in asinh(x) for mu calculation
  real x =
      0.5 / drParameters.rsSr0 * std::exp(localStateVariable / localA) * localSlipRateMagnitude;
  return localA * misc::asinh(x);
}

real FastVelocityWeakeningLaw::updateMuDerivative(unsigned int ltsFace,
                                                  unsigned int pointIndex,
                                                  real localSlipRateMagnitude,
                                                  real localStateVariable) {
  real localA = a[ltsFace][pointIndex];
  real c = 0.5 / drParameters.rsSr0 * std::exp(localStateVariable / localA);
  return localA * c / std::sqrt(misc::power<2>(localSlipRateMagnitude * c) + 1);
}

std::array<real, misc::numPaddedPoints> FastVelocityWeakeningLaw::resampleStateVar(
    std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned int ltsFace) {
  std::array<real, misc::numPaddedPoints> deltaStateVar = {0};
  std::array<real, misc::numPaddedPoints> resampledDeltaStateVar = {0};
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; ++pointIndex) {
    deltaStateVar[pointIndex] =
        stateVariableBuffer[pointIndex] - this->stateVariable[ltsFace][pointIndex];
  }
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resampleM = init::resample::Values;
  resampleKrnl.resamplePar = deltaStateVar.data();
  resampleKrnl.resampledPar = resampledDeltaStateVar.data(); // output from execute
  resampleKrnl.execute();

  return resampledDeltaStateVar;
}

void RateAndStateThermalPressurizationLaw::initializeTP(
    seissol::Interoperability& eInteroperability) {
  eInteroperability.getDynRupTP(tpGrid, tpDFinv);
}

void RateAndStateThermalPressurizationLaw::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  FastVelocityWeakeningLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);
  temperature = layerData.var(concreteLts->temperature);
  pressure = layerData.var(concreteLts->pressure);
  tpTheta = layerData.var(concreteLts->tpTheta);
  tpSigma = layerData.var(concreteLts->tpSigma);
  tpHalfWidthShearZone = layerData.var(concreteLts->tpHalfWidthShearZone);
  alphaHy = layerData.var(concreteLts->alphaHy);
}

void RateAndStateThermalPressurizationLaw::setInitialFluidPressureHook(
    std::array<real, misc::numPaddedPoints>& fluidPressure, unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    fluidPressure[pointIndex] = pressure[ltsFace][pointIndex];
  }
}

void RateAndStateThermalPressurizationLaw::calcFluidPressureHook(
    std::array<real, misc::numPaddedPoints>& fluidPressure,
    FaultStresses& faultStresses,
    bool saveTmpInTP,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

    // compute fault strength
    auto normalStress = faultStresses.normalStress[timeIndex][pointIndex] +
                        initialStressInFaultCS[ltsFace][pointIndex][0] - fluidPressure[pointIndex];
    faultStrength[pointIndex] =
        -mu[ltsFace][pointIndex] * std::min(static_cast<real>(0.0), normalStress);

    for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < numberOfTPGridPoints;
         tpGridPointIndex++) {
      //! recover original values as it gets overwritten in the ThermalPressure routine
      thetaTmp[tpGridPointIndex] = tpTheta[ltsFace][pointIndex][tpGridPointIndex];
      sigmaTmp[tpGridPointIndex] = tpSigma[ltsFace][pointIndex][tpGridPointIndex];
    }
    //! use Theta/Sigma from last call in this update, dt/2 and new SR from NS
    updateTemperatureAndPressure(pointIndex, timeIndex, ltsFace);

    fluidPressure[pointIndex] = pressure[ltsFace][pointIndex];
    if (saveTmpInTP) {
      for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < numberOfTPGridPoints;
           tpGridPointIndex++) {
        tpTheta[ltsFace][pointIndex][tpGridPointIndex] = thetaTmp[tpGridPointIndex];
        tpSigma[ltsFace][pointIndex][tpGridPointIndex] = sigmaTmp[tpGridPointIndex];
      }
    }
  }
}

void RateAndStateThermalPressurizationLaw::updateTemperatureAndPressure(unsigned int pointIndex,
                                                                        unsigned int timeIndex,
                                                                        unsigned int ltsFace) {
  real localTemperature = 0.0;
  real localPressure = 0.0;

  real tauV = faultStrength[pointIndex] *
              slipRateMagnitude[ltsFace][pointIndex]; //! fault strenght*slip rate
  real lambdaPrime = drParameters.tpLambda * drParameters.alphaTh /
                     (alphaHy[ltsFace][pointIndex] - drParameters.alphaTh);

  real tmp[numberOfTPGridPoints]{};
  real omegaTheta[numberOfTPGridPoints]{};
  real omegaSigma[numberOfTPGridPoints]{};
  real thetaCurrent[numberOfTPGridPoints]{};
  real sigmaCurrent[numberOfTPGridPoints]{};
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < numberOfTPGridPoints;
       tpGridPointIndex++) {
    //! Gaussian shear zone in spectral domain, normalized by w
    tmp[tpGridPointIndex] =
        misc::power<2>(tpGrid[tpGridPointIndex] / tpHalfWidthShearZone[ltsFace][pointIndex]);
    //! 1. Calculate diffusion of the field at previous timestep

    //! temperature
    thetaCurrent[tpGridPointIndex] =
        thetaTmp[tpGridPointIndex] *
        exp(-drParameters.alphaTh * deltaT[timeIndex] * tmp[tpGridPointIndex]);
    //! pore pressure + lambda'*temp
    sigmaCurrent[tpGridPointIndex] =
        sigmaTmp[tpGridPointIndex] *
        exp(-alphaHy[ltsFace][pointIndex] * deltaT[timeIndex] * tmp[tpGridPointIndex]);

    //! 2. Add current contribution and get new temperature
    omegaTheta[tpGridPointIndex] =
        heatSource(tmp[tpGridPointIndex], drParameters.alphaTh, tpGridPointIndex, timeIndex);
    thetaTmp[tpGridPointIndex] =
        thetaCurrent[tpGridPointIndex] + (tauV / drParameters.rhoC) * omegaTheta[tpGridPointIndex];
    omegaSigma[tpGridPointIndex] = heatSource(
        tmp[tpGridPointIndex], alphaHy[ltsFace][pointIndex], tpGridPointIndex, timeIndex);
    sigmaTmp[tpGridPointIndex] =
        sigmaCurrent[tpGridPointIndex] + ((drParameters.tpLambda + lambdaPrime) * tauV) /
                                             (drParameters.rhoC) * omegaSigma[tpGridPointIndex];

    //! 3. Recover temperature and pressure using inverse Fourier
    //! transformation with the calculated fourier coefficients

    //! new contribution
    localTemperature += (tpDFinv[tpGridPointIndex] / tpHalfWidthShearZone[ltsFace][pointIndex]) *
                        thetaTmp[tpGridPointIndex];
    localPressure += (tpDFinv[tpGridPointIndex] / tpHalfWidthShearZone[ltsFace][pointIndex]) *
                     sigmaTmp[tpGridPointIndex];
  }
  // Update pore pressure change (sigma = pore pressure + lambda'*temp)
  // In the BIEM code (Lapusta) they use T without initial value
  localPressure = localPressure - lambdaPrime * localTemperature;

  // Temp and pore pressure change at single GP on the fault + initial values
  temperature[ltsFace][pointIndex] = localTemperature + drParameters.initialTemperature;
  pressure[ltsFace][pointIndex] = -localPressure + drParameters.initialPressure;
}

real RateAndStateThermalPressurizationLaw::heatSource(real tmp,
                                                      real alpha,
                                                      unsigned int tpGridPointIndex,
                                                      unsigned int timeIndex) {
  //! original function in spatial domain
  //! omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/TP_halfWidthShearZone).^2);
  //! function in the wavenumber domain *including additional factors in front of the heat source
  //! function* omega =
  //! 1/(*alpha*Dwn**2**(sqrt(2.0*pi))*exp(-0.5*(Dwn*TP_halfWidthShearZone)**2)*(1-exp(-alpha**dt**tmp))
  //! inserting Dwn/TP_halfWidthShearZone (scaled) for Dwn cancels out TP_halfWidthShearZone
  return 1.0 / (alpha * tmp * (sqrt(2.0 * M_PI))) *
         exp(-0.5 * misc::power<2>(tpGrid[tpGridPointIndex])) *
         (1.0 - exp(-alpha * deltaT[timeIndex] * tmp));
}
} // namespace seissol::dr::friction_law

#include "FastVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law {
void FastVelocityWeakeningLaw::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                                  seissol::initializers::DynamicRupture* dynRup,
                                                  real fullUpdateTime) {
  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);

  averagedSlip = layerData.var(concreteLts->averagedSlip);
  srW = layerData.var(concreteLts->rs_srW);
}

real FastVelocityWeakeningLaw::updateStateVariable(int pointIndex,
                                                   unsigned int face,
                                                   real stateVarReference,
                                                   real timeIncrement,
                                                   real localSlipRate) {
  double fw = drParameters.mu_w;
  double localSrW = srW[face][pointIndex];
  double localA = a[face][pointIndex];
  double localSl0 = sl0[face][pointIndex];

  // low-velocity steady state friction coefficient
  real flv =
      drParameters.rs_f0 - (drParameters.rs_b - localA) * log(localSlipRate / drParameters.rs_sr0);
  // steady state friction coefficient
  real fss = fw + (flv - fw) / pow(1.0 + std::pow(localSlipRate / localSrW, 8.0), 1.0 / 8.0);
  // steady-state state variable
  // For compiling reasons we write SINH(X)=(EXP(X)-EXP(-X))/2
  real SVss = localA * log(2.0 * drParameters.rs_sr0 / localSlipRate *
                           (exp(fss / localA) - exp(-fss / localA)) / 2.0);

  // exact integration of dSV/dt DGL, assuming constant V over integration step

  real exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
  real localStateVariable = SVss * (1.0 - exp1) + exp1 * stateVarReference;
  assert(!(std::isnan(localStateVariable) && pointIndex < numberOfPoints) && "NaN detected");
  return localStateVariable;
}

real FastVelocityWeakeningLaw::updateMu(unsigned int ltsFace,
                                        unsigned int pointIndex,
                                        real localStateVariable) {
  // mu = a * arcsinh ( V / (2*V_0) * exp (psi / a))
  real localA = a[ltsFace][pointIndex];
  // x in asinh(x) for mu calculation
  real x = 0.5 / drParameters.rs_sr0 * std::exp(localStateVariable / localA) *
           slipRateMagnitude[ltsFace][pointIndex];
  // asinh(x)=log(x+sqrt(x^2+1))
  return localA * std::log(x + std::sqrt(std::pow(x, 2) + 1.0));
}

void RateAndStateThermalPressurizationLaw::initializeTP(
    seissol::Interoperability& e_interoperability) {
  e_interoperability.getDynRupTP(TP_grid, TP_DFinv);
}

void RateAndStateThermalPressurizationLaw::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  FastVelocityWeakeningLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateThermalPressurisation*>(dynRup);
  temperature = layerData.var(concreteLts->temperature);
  pressure = layerData.var(concreteLts->pressure);
  TP_Theta = layerData.var(concreteLts->TP_theta);
  TP_sigma = layerData.var(concreteLts->TP_sigma);
  TP_halfWidthShearZone = layerData.var(concreteLts->TP_halfWidthShearZone);
  alphaHy = layerData.var(concreteLts->alphaHy);
}

void RateAndStateThermalPressurizationLaw::hookSetInitialP_f(std::array<real, numPaddedPoints>& P_f,
                                                             unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    P_f[pointIndex] = pressure[ltsFace][pointIndex];
  }
}

void RateAndStateThermalPressurizationLaw::hookCalcP_f(std::array<real, numPaddedPoints>& P_f,
                                                       FaultStresses& faultStresses,
                                                       bool saveTmpInTP,
                                                       unsigned int timeIndex,
                                                       unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

    // compute fault strength (Sh)
    faultStrength[pointIndex] = -mu[ltsFace][pointIndex] *
                                (faultStresses.NormalStressGP[timeIndex][pointIndex] +
                                 initialStressInFaultCS[ltsFace][pointIndex][0] - P_f[pointIndex]);

    for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < numberOfTPGridPoints;
         tpGridPointIndex++) {
      //! recover original values as it gets overwritten in the ThermalPressure routine
      Theta_tmp[tpGridPointIndex] = TP_Theta[ltsFace][pointIndex][tpGridPointIndex];
      Sigma_tmp[tpGridPointIndex] = TP_sigma[ltsFace][pointIndex][tpGridPointIndex];
    }
    //! use Theta/Sigma from last call in this update, dt/2 and new SR from NS
    updateTemperatureAndPressure(pointIndex, timeIndex, ltsFace);

    P_f[pointIndex] = pressure[ltsFace][pointIndex];
    if (saveTmpInTP) {
      for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < numberOfTPGridPoints;
           tpGridPointIndex++) {
        TP_Theta[ltsFace][pointIndex][tpGridPointIndex] = Theta_tmp[tpGridPointIndex];
        TP_sigma[ltsFace][pointIndex][tpGridPointIndex] = Sigma_tmp[tpGridPointIndex];
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
  real lambdaPrime = drParameters.tP_lambda * drParameters.alpha_th /
                     (alphaHy[ltsFace][pointIndex] - drParameters.alpha_th);

  real tmp[numberOfTPGridPoints]{};
  real omegaTheta[numberOfTPGridPoints]{};
  real omegaSigma[numberOfTPGridPoints]{};
  real thetaCurrent[numberOfTPGridPoints]{};
  real sigmaCurrent[numberOfTPGridPoints]{};
  for (unsigned int tpGridPointIndex = 0; tpGridPointIndex < numberOfTPGridPoints;
       tpGridPointIndex++) {
    //! Gaussian shear zone in spectral domain, normalized by w
    tmp[tpGridPointIndex] =
        std::pow(TP_grid[tpGridPointIndex] / TP_halfWidthShearZone[ltsFace][pointIndex], 2);
    //! 1. Calculate diffusion of the field at previous timestep

    //! temperature
    thetaCurrent[tpGridPointIndex] =
        Theta_tmp[tpGridPointIndex] *
        exp(-drParameters.alpha_th * deltaT[timeIndex] * tmp[tpGridPointIndex]);
    //! pore pressure + lambda'*temp
    sigmaCurrent[tpGridPointIndex] =
        Sigma_tmp[tpGridPointIndex] *
        exp(-alphaHy[ltsFace][pointIndex] * deltaT[timeIndex] * tmp[tpGridPointIndex]);

    //! 2. Add current contribution and get new temperature
    omegaTheta[tpGridPointIndex] =
        heatSource(tmp[tpGridPointIndex], drParameters.alpha_th, tpGridPointIndex, timeIndex);
    Theta_tmp[tpGridPointIndex] =
        thetaCurrent[tpGridPointIndex] + (tauV / drParameters.rho_c) * omegaTheta[tpGridPointIndex];
    omegaSigma[tpGridPointIndex] = heatSource(
        tmp[tpGridPointIndex], alphaHy[ltsFace][pointIndex], tpGridPointIndex, timeIndex);
    Sigma_tmp[tpGridPointIndex] =
        sigmaCurrent[tpGridPointIndex] + ((drParameters.tP_lambda + lambdaPrime) * tauV) /
                                             (drParameters.rho_c) * omegaSigma[tpGridPointIndex];

    //! 3. Recover temperature and pressure using inverse Fourier
    //! transformation with the calculated fourier coefficients

    //! new contribution
    localTemperature += (TP_DFinv[tpGridPointIndex] / TP_halfWidthShearZone[ltsFace][pointIndex]) *
                        Theta_tmp[tpGridPointIndex];
    localPressure += (TP_DFinv[tpGridPointIndex] / TP_halfWidthShearZone[ltsFace][pointIndex]) *
                     Sigma_tmp[tpGridPointIndex];
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
         exp(-0.5 * std::pow(TP_grid[tpGridPointIndex], 2)) *
         (1.0 - exp(-alpha * deltaT[timeIndex] * tmp));
}
} // namespace seissol::dr::friction_law

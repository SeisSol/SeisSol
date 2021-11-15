#include "RateAndState.h"

namespace seissol::dr::friction_law {
void RateAndStateFastVelocityWeakeningLaw::copyLtsTreeToLocalRS(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  // maybe change later to const_cast?
  auto concreteLts =
      dynamic_cast<seissol::initializers::LTS_RateAndStateFastVelocityWeakening*>(dynRup);
  nucleationStressInFaultCS = layerData.var(concreteLts->nucleationStressInFaultCS);

  sl0 = layerData.var(concreteLts->rs_sl0);
  a = layerData.var(concreteLts->rs_a);
  srW = layerData.var(concreteLts->rs_srW);
  ds = layerData.var(concreteLts->ds);
  averagedSlip = layerData.var(concreteLts->averagedSlip);
  stateVariable = layerData.var(concreteLts->stateVariable);
  dynStressTime = layerData.var(concreteLts->dynStressTime);
}

void RateAndStateFastVelocityWeakeningLaw::preCalcTime() {
  dt = 0;
  for (int timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
    dt += deltaT[timeIndex];
  }
  if (m_fullUpdateTime <= drParameters.t0) {
    gNuc = calcSmoothStepIncrement(m_fullUpdateTime, dt);
  }
}

void RateAndStateFastVelocityWeakeningLaw::setInitialValues(
    std::array<real, numPaddedPoints>& localStateVariable, unsigned int ltsFace) {
  if (m_fullUpdateTime <= drParameters.t0) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      for (int i = 0; i < 6; i++) {
        initialStressInFaultCS[ltsFace][pointIndex][i] +=
            nucleationStressInFaultCS[ltsFace][pointIndex][i] * gNuc;
      }
    }
  } // end If-Tnuc

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    localStateVariable[pointIndex] = stateVariable[ltsFace][pointIndex];
  }
}

void RateAndStateFastVelocityWeakeningLaw::calcInitialSlipRate(
    std::array<real, numPaddedPoints>& totalShearStressYZ,
    FaultStresses& faultStresses,
    std::array<real, numPaddedPoints>& stateVarZero,
    std::array<real, numPaddedPoints>& localStateVariable,
    std::array<real, numPaddedPoints>& temporarySlipRate,
    unsigned int timeIndex,
    unsigned int ltsFace) {

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

    // friction develops as                    mu = a * arcsinh[ V/(2*V0) * exp(SV/a) ]
    // state variable SV develops as          dSV / dt = -(V - L) * (SV - SV_ss)
    //                                        SV_ss = a * ln[ 2*V0/V * sinh(mu_ss/a) ]
    //                                        mu_ss = mu_w + [mu_lv - mu_w] / [ 1 + (V/Vw)^8 ] ^
    //                                        (1/8) ] mu_lv = mu_0 - (b-a) ln (V/V0)

    totalShearStressYZ[pointIndex] =
        std::sqrt(std::pow(initialStressInFaultCS[ltsFace][pointIndex][3] +
                               faultStresses.XYStressGP[timeIndex][pointIndex],
                           2) +
                  std::pow(initialStressInFaultCS[ltsFace][pointIndex][5] +
                               faultStresses.XZStressGP[timeIndex][pointIndex],
                           2));

    // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
    // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )
    // Careful, the state variable must always be corrected using stateVar0 and not
    // localStateVariable!
    stateVarZero[pointIndex] = localStateVariable[pointIndex];

    // The following process is adapted from that described by Kaneko et al. (2008)
    slipRateMagnitude[ltsFace][pointIndex] =
        std::sqrt(std::pow(slipRateStrike[ltsFace][pointIndex], 2) +
                  std::pow(slipRateDip[ltsFace][pointIndex], 2));
    slipRateMagnitude[ltsFace][pointIndex] =
        std::max(almostZero, slipRateMagnitude[ltsFace][pointIndex]);
    temporarySlipRate[pointIndex] = slipRateMagnitude[ltsFace][pointIndex];
  } // End of pointIndex-loop
}

void RateAndStateFastVelocityWeakeningLaw::updateStateVariableIterative(
    bool& has_converged,
    std::array<real, numPaddedPoints>& stateVarZero,
    std::array<real, numPaddedPoints>& SR_tmp,
    std::array<real, numPaddedPoints>& localStateVariable,
    std::array<real, numPaddedPoints>& P_f,
    std::array<real, numPaddedPoints>& normalStress,
    std::array<real, numPaddedPoints>& TotalShearStressYZ,
    std::array<real, numPaddedPoints>& SRtest,
    FaultStresses& faultStresses,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    // fault strength using LocMu and P_f from previous timestep/iteration
    // 1.update SV using Vold from the previous time step
    updateStateVariable(pointIndex,
                        ltsFace,
                        stateVarZero[pointIndex],
                        deltaT[timeIndex],
                        SR_tmp[pointIndex],
                        localStateVariable[pointIndex]);
    normalStress[pointIndex] = faultStresses.NormalStressGP[timeIndex][pointIndex] +
                               initialStressInFaultCS[ltsFace][pointIndex][0] - P_f[pointIndex];
  } // End of pointIndex-loop

  // 2. solve for Vnew , applying the Newton-Raphson algorithm
  // effective normal stress including initial stresses and pore fluid pressure
  has_converged = IterativelyInvertSR(
      ltsFace, numberSlipRateUpdates, localStateVariable, normalStress, TotalShearStressYZ, SRtest);

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    // 3. update theta, now using V=(Vnew+Vold)/2
    // For the next SV update, use the mean slip rate between the initial guess and the one found
    // (Kaneko 2008, step 6)
    SR_tmp[pointIndex] = 0.5 * (slipRateMagnitude[ltsFace][pointIndex] + fabs(SRtest[pointIndex]));

    // 4. solve again for Vnew
    slipRateMagnitude[ltsFace][pointIndex] = fabs(SRtest[pointIndex]);

    // update LocMu
    updateMu(ltsFace, pointIndex, localStateVariable[pointIndex]);
  } // End of pointIndex-loop
}

void RateAndStateFastVelocityWeakeningLaw::executeIfNotConverged(
    std::array<real, numPaddedPoints>& localStateVariable, unsigned ltsFace) {
  real tmp = 0.5 / drParameters.rs_sr0 * exp(localStateVariable[0] / a[ltsFace][0]) *
             slipRateMagnitude[ltsFace][0];
  //! logError(*) 'nonConvergence RS Newton', time
  std::cout << "nonConvergence RS Newton, time: " << m_fullUpdateTime << std::endl;
  assert(!std::isnan(tmp) && "nonConvergence RS Newton");
}

void RateAndStateFastVelocityWeakeningLaw::calcSlipRateAndTraction(
    std::array<real, numPaddedPoints>& stateVarZero,
    std::array<real, numPaddedPoints>& SR_tmp,
    std::array<real, numPaddedPoints>& localStateVariable,
    std::array<real, numPaddedPoints>& normalStress,
    std::array<real, numPaddedPoints>& TotalShearStressYZ,
    std::array<real, numPaddedPoints>& tmpSlip,
    real deltaStateVar[numPaddedPoints],
    FaultStresses& faultStresses,
    unsigned int timeIndex,
    unsigned int ltsFace) {
  std::array<real, numPaddedPoints> LocslipRateMagnitude{0};

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //! SV from mean slip rate in tmp
    updateStateVariable(pointIndex,
                        ltsFace,
                        stateVarZero[pointIndex],
                        deltaT[timeIndex],
                        SR_tmp[pointIndex],
                        localStateVariable[pointIndex]);

    //! update LocMu for next strength determination, only needed for last update
    updateMu(ltsFace, pointIndex, localStateVariable[pointIndex]);

    //! update stress change
    tractionXY[ltsFace][pointIndex] = -((initialStressInFaultCS[ltsFace][pointIndex][3] +
                                         faultStresses.XYStressGP[timeIndex][pointIndex]) /
                                        TotalShearStressYZ[pointIndex]) *
                                      mu[ltsFace][pointIndex] * normalStress[pointIndex];
    tractionXZ[ltsFace][pointIndex] = -((initialStressInFaultCS[ltsFace][pointIndex][5] +
                                         faultStresses.XZStressGP[timeIndex][pointIndex]) /
                                        TotalShearStressYZ[pointIndex]) *
                                      mu[ltsFace][pointIndex] * normalStress[pointIndex];
    tractionXY[ltsFace][pointIndex] -= initialStressInFaultCS[ltsFace][pointIndex][3];
    tractionXZ[ltsFace][pointIndex] -= initialStressInFaultCS[ltsFace][pointIndex][5];

    // Compute slip
    //! ABS of locSlipRate removed as it would be the accumulated slip that is usually not needed in
    //! the solver, see linear slip weakening
    slip[ltsFace][pointIndex] += slipRateMagnitude[ltsFace][pointIndex] * deltaT[timeIndex];

    //! Update slip rate (notice that locSlipRate(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate
    //! caused by a free surface!)
    slipRateStrike[ltsFace][pointIndex] =
        -impAndEta[ltsFace].inv_eta_s *
        (tractionXY[ltsFace][pointIndex] - faultStresses.XYStressGP[timeIndex][pointIndex]);
    slipRateDip[ltsFace][pointIndex] =
        -impAndEta[ltsFace].inv_eta_s *
        (tractionXZ[ltsFace][pointIndex] - faultStresses.XZStressGP[timeIndex][pointIndex]);

    //! TU 07.07.16: correct slipRateStrike and slipRateDip to avoid numerical errors
    LocslipRateMagnitude[pointIndex] = sqrt(std::pow(slipRateStrike[ltsFace][pointIndex], 2) +
                                            std::pow(slipRateDip[ltsFace][pointIndex], 2));
    if (LocslipRateMagnitude[pointIndex] != 0) {
      slipRateStrike[ltsFace][pointIndex] *=
          slipRateMagnitude[ltsFace][pointIndex] / LocslipRateMagnitude[pointIndex];
      slipRateDip[ltsFace][pointIndex] *=
          slipRateMagnitude[ltsFace][pointIndex] / LocslipRateMagnitude[pointIndex];
    }
    tmpSlip[pointIndex] += LocslipRateMagnitude[pointIndex] * deltaT[timeIndex];

    slipStrike[ltsFace][pointIndex] += slipRateStrike[ltsFace][pointIndex] * deltaT[timeIndex];
    slipDip[ltsFace][pointIndex] += slipRateDip[ltsFace][pointIndex] * deltaT[timeIndex];

    //! Save traction for flux computation
    faultStresses.XYTractionResultGP[timeIndex][pointIndex] = tractionXY[ltsFace][pointIndex];
    faultStresses.XZTractionResultGP[timeIndex][pointIndex] = tractionXZ[ltsFace][pointIndex];

    // Could be outside TimeLoop, since only last time result is used later
    deltaStateVar[pointIndex] = localStateVariable[pointIndex] - stateVariable[ltsFace][pointIndex];
  } // End of BndGP-loop
}

void RateAndStateFastVelocityWeakeningLaw::resampleStateVar(real deltaStateVar[numPaddedPoints],
                                                            unsigned int ltsFace) {
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resampleM = init::resample::Values;
  real resampledDeltaStateVariable[numPaddedPoints];
  resampleKrnl.resamplePar = deltaStateVar;
  resampleKrnl.resampledPar = resampledDeltaStateVariable; // output from execute
  resampleKrnl.execute();

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    // write back State Variable to lts tree
    stateVariable[ltsFace][pointIndex] =
        stateVariable[ltsFace][pointIndex] + resampledDeltaStateVariable[pointIndex];
  }
}

void RateAndStateFastVelocityWeakeningLaw::saveDynamicStressOutput(unsigned int face) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

    if (ruptureTime[face][pointIndex] > 0.0 && ruptureTime[face][pointIndex] <= m_fullUpdateTime &&
        ds[pointIndex] &&
        mu[face][pointIndex] <=
            (drParameters.mu_w + 0.05 * (drParameters.rs_f0 - drParameters.mu_w))) {
      dynStressTime[face][pointIndex] = m_fullUpdateTime;
      ds[face][pointIndex] = false;
    }
  }
}

void RateAndStateFastVelocityWeakeningLaw::hookSetInitialP_f(std::array<real, numPaddedPoints>& P_f,
                                                             unsigned int ltsFace) {
  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    P_f[pointIndex] = 0.0;
  }
}

void RateAndStateFastVelocityWeakeningLaw::hookCalcP_f(std::array<real, numPaddedPoints>& P_f,
                                                       FaultStresses& faultStresses,
                                                       bool saveTmpInTP,
                                                       unsigned int timeIndex,
                                                       unsigned int ltsFace) {}

void RateAndStateFastVelocityWeakeningLaw::updateStateVariable(int pointIndex,
                                                               unsigned int face,
                                                               real SV0,
                                                               real time_inc,
                                                               real& SR_tmp,
                                                               real& localStateVariable) {
  double flv, fss, SVss;
  double fw = drParameters.mu_w;
  double localSrW = srW[face][pointIndex];
  double localA = a[face][pointIndex];
  double localSl0 = sl0[face][pointIndex];
  double exp1;

  // low-velocity steady state friction coefficient
  flv = drParameters.rs_f0 - (drParameters.rs_b - localA) * log(SR_tmp / drParameters.rs_sr0);
  // steady state friction coefficient
  fss = fw + (flv - fw) / pow(1.0 + std::pow(SR_tmp / localSrW, 8.0), 1.0 / 8.0);
  // steady-state state variable
  // For compiling reasons we write SINH(X)=(EXP(X)-EXP(-X))/2
  SVss = localA *
         log(2.0 * drParameters.rs_sr0 / SR_tmp * (exp(fss / localA) - exp(-fss / localA)) / 2.0);

  // exact integration of dSV/dt DGL, assuming constant V over integration step

  exp1 = exp(-SR_tmp * (time_inc / localSl0));
  localStateVariable = SVss * (1.0 - exp1) + exp1 * SV0;

  assert(!(std::isnan(localStateVariable) && pointIndex < numberOfPoints) && "NaN detected");
}

bool RateAndStateFastVelocityWeakeningLaw::IterativelyInvertSR(
    unsigned int ltsFace,
    int nSRupdates,
    std::array<real, numPaddedPoints>& localStateVariable,
    std::array<real, numPaddedPoints>& normalStress,
    std::array<real, numPaddedPoints>& shearStress,
    std::array<real, numPaddedPoints>& slipRateTest) {

  real tmp[numPaddedPoints], tmp2[numPaddedPoints], tmp3[numPaddedPoints], mu_f[numPaddedPoints],
      dmu_f[numPaddedPoints], NR[numPaddedPoints], dNR[numPaddedPoints];
  // double AlmostZero = 1e-45;
  bool has_converged = false;

  //! solve for Vnew = SR , applying the Newton-Raphson algorithm
  //! SR fulfills g(SR)=f(SR)
  //!-> find root of NR=f-g using a Newton-Raphson algorithm with dNR = d(NR)/d(SR)
  //! SR_{i+1}=SR_i-( NR_i / dNR_i )
  //!
  //!        equalize:
  //!         g = SR*MU/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
  //!         f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
  //!  where mu = friction coefficient, dependening on the RSF law used

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    //! first guess = SR value of the previous step
    slipRateTest[pointIndex] = slipRateMagnitude[ltsFace][pointIndex];
    tmp[pointIndex] =
        0.5 / drParameters.rs_sr0 * exp(localStateVariable[pointIndex] / a[ltsFace][pointIndex]);
  }

  for (int i = 0; i < nSRupdates; i++) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

      //! f = ( tmp2 * ABS(LocP+P_0)- ABS(S_0))*(S_0)/ABS(S_0)
      //! g = slipRateTest * 1.0/(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) + ABS(ShTest)
      //! for compiling reasons ASINH(X)=LOG(X+SQRT(X^2+1))

      //! calculate friction coefficient
      tmp2[pointIndex] = tmp[pointIndex] * slipRateTest[pointIndex];
      mu_f[pointIndex] = a[ltsFace][pointIndex] *
                         log(tmp2[pointIndex] + sqrt(std::pow(tmp2[pointIndex], 2) + 1.0));
      dmu_f[pointIndex] =
          a[ltsFace][pointIndex] / sqrt(1.0 + std::pow(tmp2[pointIndex], 2)) * tmp[pointIndex];
      NR[pointIndex] =
          -impAndEta[ltsFace].inv_eta_s *
              (fabs(normalStress[pointIndex]) * mu_f[pointIndex] - shearStress[pointIndex]) -
          slipRateTest[pointIndex];
    }

    has_converged = true;

    // max element of NR must be smaller then aTolF
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      if (fabs(NR[pointIndex]) >= aTolF) {
        has_converged = false;
        break;
      }
    }
    if (has_converged) {
      return has_converged;
    }
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

      //! derivative of NR
      dNR[pointIndex] =
          -impAndEta[ltsFace].inv_eta_s * (fabs(normalStress[pointIndex]) * dmu_f[pointIndex]) -
          1.0;
      //! ratio
      tmp3[pointIndex] = NR[pointIndex] / dNR[pointIndex];

      //! update slipRateTest
      slipRateTest[pointIndex] = std::max(almostZero, slipRateTest[pointIndex] - tmp3[pointIndex]);
    }
  }
  return has_converged;
}

bool RateAndStateFastVelocityWeakeningLaw::IterativelyInvertSR_Brent(
    unsigned int ltsFace,
    int nSRupdates,
    std::array<real, numPaddedPoints>& localStateVariable,
    std::array<real, numPaddedPoints>& n_stress,
    std::array<real, numPaddedPoints>& sh_stress,
    std::array<real, numPaddedPoints>& SRtest) {
  std::function<double(double, int)> F;
  double tol = 1e-30;

  real* localA = a[ltsFace];
  double rs_sr0_ = drParameters.rs_sr0;
  double invEta = impAndEta[ltsFace].inv_eta_s;

  F = [invEta, &sh_stress, n_stress, localA, localStateVariable, rs_sr0_](double SR,
                                                                          int pointIndex) {
    double tmp = 0.5 / rs_sr0_ * std::exp(localStateVariable[pointIndex] / localA[pointIndex]) * SR;
    double mu_f = localA[pointIndex] * std::log(tmp + std::sqrt(std::pow(tmp, 2) + 1.0));
    return -invEta * (fabs(n_stress[pointIndex]) * mu_f - sh_stress[pointIndex]) - SR;
  };

  for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
    // TODO: use better boundaries?
    double a = slipRateMagnitude[ltsFace][pointIndex] -
               impAndEta[ltsFace].inv_eta_s * sh_stress[pointIndex];
    double b = slipRateMagnitude[ltsFace][pointIndex] +
               impAndEta[ltsFace].inv_eta_s * sh_stress[pointIndex];

    double eps = std::numeric_limits<double>::epsilon();
    double Fa = F(a, pointIndex);
    // if(std::isinf(Fa)){
    //  Fa = std::numeric_limits<double>::max();
    //}
    double Fb = F(b, pointIndex);
    assert(std::copysign(Fa, Fb) != Fa); // Fa and Fb have different signs
    double c = a;
    double Fc = Fa;
    double d = b - a;
    double e = d;
    while (Fb != 0.0) {
      if (std::copysign(Fb, Fc) == Fb) {
        c = a;
        Fc = Fa;
        d = b - a;
        e = d;
      }
      if (std::fabs(Fc) < std::fabs(Fb)) {
        a = b;
        b = c;
        c = a;
        Fa = Fb;
        Fb = Fc;
        Fc = Fa;
      }
      // Convergence test
      double xm = 0.5 * (c - b);
      double tol1 = 2.0 * eps * std::fabs(b) + 0.5 * tol;
      if (std::fabs(xm) <= tol1 || Fb == 0.0) {
        break;
      }
      if (std::fabs(e) < tol1 || std::fabs(Fa) <= std::fabs(Fb)) {
        // bisection
        d = xm;
        e = d;
      } else {
        double s = Fb / Fa;
        double p, q;
        if (a != c) {
          // linear interpolation
          q = Fa / Fc;
          double r = Fb / Fc;
          p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
          q = (q - 1.0) * (r - 1.0) * (s - 1.0);
        } else {
          // inverse quadratic interpolation
          p = 2.0 * xm * s;
          q = 1.0 - s;
        }
        if (p > 0) {
          q = -q;
        } else {
          p = -p;
        }
        if (2.0 * p < 3.0 * xm * q - std::fabs(tol1 * q) && p < std::fabs(0.5 * e * q)) {
          e = d;
          d = p / q;
        } else {
          // bisection
          d = xm;
          e = d;
        }
      }
      a = b;
      Fa = Fb;
      if (std::fabs(d) > tol1) {
        b += d;
      } else {
        b += std::copysign(tol1, xm);
      }
      Fb = F(b, pointIndex);
    }
    SRtest[pointIndex] = b;
  }
  return true;
}
void RateAndStateFastVelocityWeakeningLaw::updateMu(unsigned int ltsFace,
                                                    unsigned int pointIndex,
                                                    real localStateVariable) {
  //! X in Asinh(x) for mu calculation
  real tmp = 0.5 / drParameters.rs_sr0 * exp(localStateVariable / a[ltsFace][pointIndex]) *
             slipRateMagnitude[ltsFace][pointIndex];
  //! mu from locSlipRate with SINH(X)=LOG(X+SQRT(X^2+1))
  mu[ltsFace][pointIndex] =
      a[ltsFace][pointIndex] * std::log(tmp + std::sqrt(std::pow(tmp, 2) + 1.0));
}

void RateAndStateThermalPressurizationLaw::initializeTP(
    seissol::Interoperability& e_interoperability) {
  e_interoperability.getDynRupTP(TP_grid, TP_DFinv);
}

void RateAndStateThermalPressurizationLaw::copyLtsTreeToLocalRS(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture* dynRup,
    real fullUpdateTime) {
  // first copy all Variables from the Base Lts dynRup tree
  RateAndStateFastVelocityWeakeningLaw::copyLtsTreeToLocalRS(layerData, dynRup, fullUpdateTime);

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

    for (unsigned int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++) {
      //! recover original values as it gets overwritten in the ThermalPressure routine
      Theta_tmp[iTP_grid_nz] = TP_Theta[ltsFace][pointIndex][iTP_grid_nz];
      Sigma_tmp[iTP_grid_nz] = TP_sigma[ltsFace][pointIndex][iTP_grid_nz];
    }
    //! use Theta/Sigma from last call in this update, dt/2 and new SR from NS
    Calc_ThermalPressure(pointIndex, timeIndex, ltsFace);

    P_f[pointIndex] = pressure[ltsFace][pointIndex];
    if (saveTmpInTP) {
      for (unsigned int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++) {
        TP_Theta[ltsFace][pointIndex][iTP_grid_nz] = Theta_tmp[iTP_grid_nz];
        TP_sigma[ltsFace][pointIndex][iTP_grid_nz] = Sigma_tmp[iTP_grid_nz];
      }
    }
  }
}

void RateAndStateThermalPressurizationLaw::Calc_ThermalPressure(unsigned int pointIndex,
                                                                unsigned int timeIndex,
                                                                unsigned int ltsFace) {
  real localTemperature = 0.0;
  real localPressure = 0.0;

  real tauV = faultStrength[pointIndex] *
              slipRateMagnitude[ltsFace][pointIndex]; //! fault strenght*slip rate
  real lambdaPrime = drParameters.tP_lambda * drParameters.alpha_th /
                     (alphaHy[ltsFace][pointIndex] - drParameters.alpha_th);

  real tmp[TP_grid_nz];
  real omegaTheta[TP_grid_nz]{};
  real omegaSigma[TP_grid_nz]{};
  real thetaCurrent[TP_grid_nz]{};
  real sigmaCurrent[TP_grid_nz]{};
  for (unsigned int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++) {
    //! Gaussian shear zone in spectral domain, normalized by w
    tmp[iTP_grid_nz] =
        std::pow(TP_grid[iTP_grid_nz] / TP_halfWidthShearZone[ltsFace][pointIndex], 2);
    //! 1. Calculate diffusion of the field at previous timestep

    //! temperature
    thetaCurrent[iTP_grid_nz] =
        Theta_tmp[iTP_grid_nz] * exp(-drParameters.alpha_th * deltaT[timeIndex] * tmp[iTP_grid_nz]);
    //! pore pressure + lambda'*temp
    sigmaCurrent[iTP_grid_nz] = Sigma_tmp[iTP_grid_nz] * exp(-alphaHy[ltsFace][pointIndex] *
                                                             deltaT[timeIndex] * tmp[iTP_grid_nz]);

    //! 2. Add current contribution and get new temperature
    omegaTheta[iTP_grid_nz] =
        heat_source(tmp[iTP_grid_nz], drParameters.alpha_th, iTP_grid_nz, timeIndex);
    Theta_tmp[iTP_grid_nz] =
        thetaCurrent[iTP_grid_nz] + (tauV / drParameters.rho_c) * omegaTheta[iTP_grid_nz];
    omegaSigma[iTP_grid_nz] =
        heat_source(tmp[iTP_grid_nz], alphaHy[ltsFace][pointIndex], iTP_grid_nz, timeIndex);
    Sigma_tmp[iTP_grid_nz] =
        sigmaCurrent[iTP_grid_nz] + ((drParameters.tP_lambda + lambdaPrime) * tauV) /
                                        (drParameters.rho_c) * omegaSigma[iTP_grid_nz];

    //! 3. Recover temperature and pressure using inverse Fourier
    //! transformation with the calculated fourier coefficients

    //! new contribution
    localTemperature += (TP_DFinv[iTP_grid_nz] / TP_halfWidthShearZone[ltsFace][pointIndex]) *
                        Theta_tmp[iTP_grid_nz];
    localPressure += (TP_DFinv[iTP_grid_nz] / TP_halfWidthShearZone[ltsFace][pointIndex]) *
                     Sigma_tmp[iTP_grid_nz];
  }
  // Update pore pressure change (sigma = pore pressure + lambda'*temp)
  // In the BIEM code (Lapusta) they use T without initial value
  localPressure = localPressure - lambdaPrime * localTemperature;

  // Temp and pore pressure change at single GP on the fault + initial values
  temperature[ltsFace][pointIndex] = localTemperature + drParameters.initialTemperature;
  pressure[ltsFace][pointIndex] = -localPressure + drParameters.iniPressure;
}

real RateAndStateThermalPressurizationLaw::heat_source(real tmp,
                                                       real alpha,
                                                       unsigned int iTP_grid_nz,
                                                       unsigned int timeIndex) {
  //! original function in spatial domain
  //! omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/TP_halfWidthShearZone).^2);
  //! function in the wavenumber domain *including additional factors in front of the heat source
  //! function* omega =
  //! 1/(*alpha*Dwn**2**(sqrt(2.0*pi))*exp(-0.5*(Dwn*TP_halfWidthShearZone)**2)*(1-exp(-alpha**dt**tmp))
  //! inserting Dwn/TP_halfWidthShearZone (scaled) for Dwn cancels out TP_halfWidthShearZone
  return 1.0 / (alpha * tmp * (sqrt(2.0 * M_PI))) * exp(-0.5 * std::pow(TP_grid[iTP_grid_nz], 2)) *
         (1.0 - exp(-alpha * deltaT[timeIndex] * tmp));
}
} // namespace seissol::dr::friction_law
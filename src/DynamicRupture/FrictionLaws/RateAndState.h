#ifndef SEISSOL_RATEANDSTATE_H
#define SEISSOL_RATEANDSTATE_H

#include "BaseFrictionLaw.h"
#include "Solver/Interoperability.h"

namespace seissol::dr::friction_law {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived>
class RateAndStateBase : public BaseFrictionLaw<RateAndStateBase<Derived>> {
  public:
  // Attributes
  // CS = coordinate system
  real (*nucleationStressInFaultCS)[numPaddedPoints][6];
  real dt = 0;
  real gNuc = 0;

  real (*a)[numPaddedPoints];
  real (*sl0)[numPaddedPoints];
  real (*stateVariable)[numPaddedPoints];

  // TU 7.07.16: if the SR is too close to zero, we will have problems (NaN)
  // as a consequence, the SR is affected the AlmostZero value when too small
  const real almostZero = 1e-45;

  /**
   * Parameters of the optimisation loops
   * absolute tolerance on the function to be optimized
   * This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be
   * the most adapted Number of iteration in the loops
   */
  const unsigned int numberSlipRateUpdates = 60;
  const unsigned int numberStateVariableUpdates = 2;
  const double newtonTolerance = 1e-8;

  using BaseFrictionLaw<RateAndStateBase<Derived>>::BaseFrictionLaw;

  protected:
  public:
  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             std::array<real, numPaddedPoints>& stateVariableBuffer,
                             std::array<real, numPaddedPoints>& strengthBuffer,
                             unsigned& ltsFace,
                             unsigned& timeIndex) {
    bool has_converged = false;

    // TODO: only allocate fluidPressure if we actually need it, i.e. when thermal pressurisation is
    // turned on.
    std::array<real, numPaddedPoints> fluidPressure{0};
    // compute initial thermal pressure (ony for FL103 TP)
    static_cast<Derived*>(this)->setInitialFluidPressureHook(fluidPressure, ltsFace);
    // compute initial slip rate and reference values
    std::array<real, numPaddedPoints> localSlipRate{0};
    std::array<real, numPaddedPoints> stateVarReference{0};
    std::array<real, numPaddedPoints> normalStress{0};
    std::array<real, numPaddedPoints> absoluteShearStress{0};
    std::tie(absoluteShearStress, normalStress, stateVarReference, localSlipRate) =
        static_cast<Derived*>(this)->calcInitialVariables(
            faultStresses, stateVariableBuffer, fluidPressure, timeIndex, ltsFace);

    // compute pressure from thermal pressurization (only FL103 TP)
    // Todo: Enable TP again
    // static_cast<Derived*>(this)->hookCalcP_f(fluidPressure, faultStresses, false, timeIndex,
    // ltsFace); compute slip rates by solving non-linear system of equations (with newton)
    this->updateStateVariableIterative(has_converged,
                                       stateVarReference,
                                       localSlipRate,
                                       stateVariableBuffer,
                                       fluidPressure,
                                       normalStress,
                                       absoluteShearStress,
                                       faultStresses,
                                       timeIndex,
                                       ltsFace);

    // check for convergence
    if (!has_converged) {
      this->executeIfNotConverged(stateVariableBuffer, ltsFace);
    }
    // compute final thermal pressure for FL103TP
    // Todo: Enable TP again
    // static_cast<Derived*>(this)->hookCalcP_f(fluidPressure, faultStresses, true, timeIndex,
    // ltsFace); compute final slip rates and traction from median value of the iterative solution
    // and initial guess
    this->calcSlipRateAndTraction(stateVarReference,
                                  localSlipRate,
                                  stateVariableBuffer,
                                  normalStress,
                                  absoluteShearStress,
                                  faultStresses,
                                  timeIndex,
                                  ltsFace);
  }

  void preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {
    // compute time increments (Gnuc)
    this->preCalcTime();

    // compute initial stress
    static_cast<Derived*>(this)->adjustInitialStress(ltsFace);
    // copy state variable from last time step
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      stateVariableBuffer[pointIndex] = stateVariable[ltsFace][pointIndex];
    }
  }

  void postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {
    std::array<real, numPaddedPoints> deltaStateVar = {0};
    for (int pointIndex = 0; pointIndex < numPaddedPoints; ++pointIndex) {
      deltaStateVar[pointIndex] =
          stateVariableBuffer[pointIndex] - this->stateVariable[ltsFace][pointIndex];
    }
    // resample state variables
    // TODO: Not for slip and ageing law
    // https://github.com/SeisSol/SeisSol/blob/master/src/Physics/Evaluate_friction_law.f90#L1085
    static_cast<Derived*>(this)->resampleStateVar(deltaStateVar, ltsFace);
  }

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
    a = layerData.var(concreteLts->rs_a);
    nucleationStressInFaultCS = layerData.var(concreteLts->nucleationStressInFaultCS);
    sl0 = layerData.var(concreteLts->rs_sl0);
    stateVariable = layerData.var(concreteLts->stateVariable);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }
  /*
   * compute time increments (Gnuc)
   */
  void preCalcTime() {
    dt = 0;
    for (int timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
      dt += this->deltaT[timeIndex];
    }
    gNuc = this->calcSmoothStepIncrement(this->m_fullUpdateTime, dt);
  }

  /*
   * compute initial stress from nucleation stress
   */
  void adjustInitialStress(unsigned int ltsFace) {
    if (this->m_fullUpdateTime <= this->drParameters.t0) {
      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
        for (int i = 0; i < 6; i++) {
          this->initialStressInFaultCS[ltsFace][pointIndex][i] +=
              nucleationStressInFaultCS[ltsFace][pointIndex][i] * gNuc;
        }
      }
    }
  }

  /*
   * Compute shear stress magnitude, effective normal stress, set reference state variable
   * (StateVarZero) Compute slip rate magnitude and set it as initial guess for iterations
   */
  std::tuple<std::array<real, numPaddedPoints>,
             std::array<real, numPaddedPoints>,
             std::array<real, numPaddedPoints>,
             std::array<real, numPaddedPoints>>
      calcInitialVariables(FaultStresses const& faultStresses,
                           std::array<real, numPaddedPoints> const& localStateVariable,
                           std::array<real, numPaddedPoints> const& fluidPressure,
                           unsigned int timeIndex,
                           unsigned int ltsFace) {
    // Careful, the state variable must always be corrected using stateVarZero and not
    // localStateVariable!
    std::array<real, numPaddedPoints> stateVarZero;
    std::copy(localStateVariable.begin(), localStateVariable.end(), stateVarZero.begin());

    std::array<real, numPaddedPoints> absoluteShearStress;
    std::array<real, numPaddedPoints> normalStress;
    std::array<real, numPaddedPoints> temporarySlipRate;

    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.XYStressGP[timeIndex][pointIndex];
      real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.XZStressGP[timeIndex][pointIndex];
      absoluteShearStress[pointIndex] = misc::magnitude(totalStressXY, totalStressXZ);

      normalStress[pointIndex] = std::min(static_cast<real>(0.0),
                                          faultStresses.NormalStressGP[timeIndex][pointIndex] +
                                              this->initialStressInFaultCS[ltsFace][pointIndex][0] -
                                              fluidPressure[pointIndex]);

      // The following process is adapted from that described by Kaneko et al. (2008)
      this->slipRateMagnitude[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(almostZero, this->slipRateMagnitude[ltsFace][pointIndex]);
      temporarySlipRate[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
    } // End of pointIndex-loop
    return {absoluteShearStress, normalStress, stateVarZero, temporarySlipRate};
  }

  void updateStateVariableIterative(bool& has_converged,
                                    std::array<real, numPaddedPoints> const& stateVarReference,
                                    std::array<real, numPaddedPoints>& localSlipRate,
                                    std::array<real, numPaddedPoints>& localStateVariable,
                                    std::array<real, numPaddedPoints> const& fluidPressure,
                                    std::array<real, numPaddedPoints> const& normalStress,
                                    std::array<real, numPaddedPoints> const& absoluteShearStress,
                                    FaultStresses const& faultStresses,
                                    unsigned int timeIndex,
                                    unsigned int ltsFace) {
    std::array<real, numPaddedPoints> testSlipRate{0};
    for (unsigned int j = 0; j < numberStateVariableUpdates; j++) {
      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
        // fault strength using friction coefficient and fluid pressure from previous
        // timestep/iteration update state variable using sliprate from the previous time step
        localStateVariable[pointIndex] =
            static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                             ltsFace,
                                                             stateVarReference[pointIndex],
                                                             this->deltaT[timeIndex],
                                                             localSlipRate[pointIndex]);
      }

      // solve for new slip rate, applying the Newton-Raphson algorithm
      // effective normal stress including initial stresses and pore fluid pressure
      has_converged = this->IterativelyInvertSR(ltsFace,
                                                numberSlipRateUpdates,
                                                localStateVariable,
                                                normalStress,
                                                absoluteShearStress,
                                                testSlipRate);

      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
        // update local slip rate, now using V=(Vnew+Vold)/2
        // For the next SV update, use the mean slip rate between the initial guess and the one
        // found (Kaneko 2008, step 6)
        localSlipRate[pointIndex] =
            0.5 * (this->slipRateMagnitude[ltsFace][pointIndex] + fabs(testSlipRate[pointIndex]));

        // 4. solve again for Vnew
        this->slipRateMagnitude[ltsFace][pointIndex] = fabs(testSlipRate[pointIndex]);

        // update friction coefficient based on new state variable and slip rate
        this->mu[ltsFace][pointIndex] = static_cast<Derived*>(this)->updateMu(
            ltsFace, pointIndex, localStateVariable[pointIndex]);
      } // End of pointIndex-loop
    }
  }

  void executeIfNotConverged(std::array<real, numPaddedPoints> const& localStateVariable,
                             unsigned ltsFace) {
    // TODO: Make this work for slow velocity weakening
    // https://github.com/SeisSol/SeisSol/blob/master/src/Physics/Evaluate_friction_law.f90#L1013
    real tmp = 0.5 / this->drParameters.rs_sr0 * exp(localStateVariable[0] / a[ltsFace][0]) *
               this->slipRateMagnitude[ltsFace][0];
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }

  void calcSlipRateAndTraction(std::array<real, numPaddedPoints> const& stateVarReference,
                               std::array<real, numPaddedPoints> const& localSlipRate,
                               std::array<real, numPaddedPoints>& localStateVariable,
                               std::array<real, numPaddedPoints> const& normalStress,
                               std::array<real, numPaddedPoints> const& absoluteShearStress,
                               FaultStresses& faultStresses,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      // SV from mean slip rate in tmp
      localStateVariable[pointIndex] =
          static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                           ltsFace,
                                                           stateVarReference[pointIndex],
                                                           this->deltaT[timeIndex],
                                                           localSlipRate[pointIndex]);

      // update LocMu for next strength determination, only needed for last update
      this->mu[ltsFace][pointIndex] = static_cast<Derived*>(this)->updateMu(
          ltsFace, pointIndex, localStateVariable[pointIndex]);
      real strength = -this->mu[ltsFace][pointIndex] * normalStress[pointIndex];
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.XYStressGP[timeIndex][pointIndex];
      real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.XZStressGP[timeIndex][pointIndex];
      // update stress change
      this->tractionXY[ltsFace][pointIndex] =
          (totalStressXY / absoluteShearStress[pointIndex]) * strength -
          this->initialStressInFaultCS[ltsFace][pointIndex][3];
      this->tractionXZ[ltsFace][pointIndex] =
          (totalStressXZ / absoluteShearStress[pointIndex]) * strength -
          this->initialStressInFaultCS[ltsFace][pointIndex][5];

      // Compute slip
      // ABS of locSlipRate removed as it would be the accumulated slip that is usually not needed
      // in the solver, see linear slip weakening
      this->slipMagnitude[ltsFace][pointIndex] +=
          this->slipRateMagnitude[ltsFace][pointIndex] * this->deltaT[timeIndex];

      // Update slip rate (notice that locSlipRate(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate
      // caused by a free surface!)
      this->slipRate1[ltsFace][pointIndex] =
          -this->impAndEta[ltsFace].inv_eta_s *
          (this->tractionXY[ltsFace][pointIndex] - faultStresses.XYStressGP[timeIndex][pointIndex]);
      this->slipRate2[ltsFace][pointIndex] =
          -this->impAndEta[ltsFace].inv_eta_s *
          (this->tractionXZ[ltsFace][pointIndex] - faultStresses.XZStressGP[timeIndex][pointIndex]);

      // TU 07.07.16: correct slipRate1 and slipRate2 to avoid numerical errors
      real locSlipRateMagnitude = misc::magnitude(this->slipRate1[ltsFace][pointIndex],
                                                  this->slipRate2[ltsFace][pointIndex]);
      if (locSlipRateMagnitude != 0) {
        this->slipRate1[ltsFace][pointIndex] *=
            this->slipRateMagnitude[ltsFace][pointIndex] / locSlipRateMagnitude;
        this->slipRate2[ltsFace][pointIndex] *=
            this->slipRateMagnitude[ltsFace][pointIndex] / locSlipRateMagnitude;
      }

      // Save traction for flux computation
      faultStresses.XYTractionResultGP[timeIndex][pointIndex] =
          this->tractionXY[ltsFace][pointIndex];
      faultStresses.XZTractionResultGP[timeIndex][pointIndex] =
          this->tractionXZ[ltsFace][pointIndex];

      // update directional slip
      this->slip1[ltsFace][pointIndex] +=
          this->slipRate1[ltsFace][pointIndex] * this->deltaT[timeIndex];
      this->slip2[ltsFace][pointIndex] +=
          this->slipRate2[ltsFace][pointIndex] * this->deltaT[timeIndex];
    } // End of BndGP-loop
  }

  void resampleStateVar(std::array<real, numPaddedPoints>& deltaStateVar, unsigned int ltsFace) {
    dynamicRupture::kernel::resampleParameter resampleKrnl;
    resampleKrnl.resampleM = init::resample::Values;
    real resampledDeltaStateVariable[numPaddedPoints];
    resampleKrnl.resamplePar = deltaStateVar.data();
    resampleKrnl.resampledPar = resampledDeltaStateVariable; // output from execute
    resampleKrnl.execute();

    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      // write back State Variable to lts tree
      stateVariable[ltsFace][pointIndex] += resampledDeltaStateVariable[pointIndex];
    }
  }

  void saveDynamicStressOutput(unsigned int face) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

      if (this->ruptureTime[face][pointIndex] > 0.0 &&
          this->ruptureTime[face][pointIndex] <= this->m_fullUpdateTime &&
          this->dynStressTimePending[pointIndex] &&
          this->mu[face][pointIndex] <=
              (this->drParameters.mu_w +
               0.05 * (this->drParameters.rs_f0 - this->drParameters.mu_w))) {
        this->dynStressTime[face][pointIndex] = this->m_fullUpdateTime;
        this->dynStressTimePending[face][pointIndex] = false;
      }
    }
  }

  void setInitialFluidPressureHook(std::array<real, numPaddedPoints>& P_f, unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      P_f[pointIndex] = 0.0;
    }
  }

  void calcFluidPressureHook(std::array<real, numPaddedPoints>& P_f,
                             FaultStresses& faultStresses,
                             bool saveTmpInTP,
                             unsigned int timeIndex,
                             unsigned int ltsFace) {}

  bool IterativelyInvertSR(unsigned int ltsFace,
                           int nSRupdates,
                           std::array<real, numPaddedPoints> const& localStateVariable,
                           std::array<real, numPaddedPoints> const& normalStress,
                           std::array<real, numPaddedPoints> const& absoluteShearStress,
                           std::array<real, numPaddedPoints>& slipRateTest) {

    real tmp[numPaddedPoints], tmp2[numPaddedPoints], tmp3[numPaddedPoints];
    real mu_f[numPaddedPoints], dmu_f[numPaddedPoints];
    real NR[numPaddedPoints], dNR[numPaddedPoints];
    // double AlmostZero = 1e-45;
    bool has_converged = false;

    // solve for Vnew = SR , applying the Newton-Raphson algorithm
    // SR fulfills g(SR)=f(SR)
    // -> find root of NR=f-g using a Newton-Raphson algorithm with dNR = d(NR)/d(SR)
    // SR_{i+1}=SR_i-( NR_i / dNR_i )
    //
    //        equalize:
    //         g = SR*MU/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
    //         f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
    //  where mu = friction coefficient, dependening on the RSF law used

    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      // first guess = SR value of the previous step
      slipRateTest[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
      // TODO tmp is different for slow velocity weakening.
      // https://github.com/SeisSol/SeisSol/blame/master/src/Physics/Evaluate_friction_law.f90#L1190
      tmp[pointIndex] = 0.5 / this->drParameters.rs_sr0 *
                        exp(localStateVariable[pointIndex] / a[ltsFace][pointIndex]);
    }

    for (int i = 0; i < nSRupdates; i++) {
      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

        // f = ( tmp2 * ABS(LocP+P_0)- ABS(S_0))*(S_0)/ABS(S_0)
        // g = slipRateTest * 1.0/(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) + ABS(ShTest)
        // for compiling reasons ASINH(X)=LOG(X+SQRT(X^2+1))

        // calculate friction coefficient
        tmp2[pointIndex] = tmp[pointIndex] * slipRateTest[pointIndex];
        mu_f[pointIndex] = a[ltsFace][pointIndex] * misc::asinh(tmp2[pointIndex]);
        dmu_f[pointIndex] =
            a[ltsFace][pointIndex] / sqrt(1.0 + misc::power<2>(tmp2[pointIndex])) * tmp[pointIndex];
        NR[pointIndex] = -this->impAndEta[ltsFace].inv_eta_s *
                             (fabs(normalStress[pointIndex]) * mu_f[pointIndex] -
                              absoluteShearStress[pointIndex]) -
                         slipRateTest[pointIndex];
      }

      has_converged = true;

      // max element of NR must be smaller than newtonTolerance
      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
        if (fabs(NR[pointIndex]) >= newtonTolerance) {
          has_converged = false;
          break;
        }
      }
      if (has_converged) {
        return has_converged;
      }
      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {

        // derivative of NR
        dNR[pointIndex] = -this->impAndEta[ltsFace].inv_eta_s *
                              (fabs(normalStress[pointIndex]) * dmu_f[pointIndex]) -
                          1.0;
        // ratio
        tmp3[pointIndex] = NR[pointIndex] / dNR[pointIndex];

        // update slipRateTest
        slipRateTest[pointIndex] =
            std::max(almostZero, slipRateTest[pointIndex] - tmp3[pointIndex]);
      }
    }
    return has_converged;
  }
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_RATEANDSTATE_H
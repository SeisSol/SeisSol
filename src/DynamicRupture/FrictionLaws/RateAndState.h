#ifndef SEISSOL_RATEANDSTATE_H
#define SEISSOL_RATEANDSTATE_H

#include "BaseFrictionLaw.h"
#include "Solver/Interoperability.h"

namespace seissol::dr::friction_law {
/**
 * General implementation of a rate and state solver
 * Methods are inherited via CRTP and must be implemented in the child class.
 */
template <class Derived, class TPMethod>
class RateAndStateBase : public BaseFrictionLaw<RateAndStateBase<Derived, TPMethod>> {
  public:
  // Attributes
  // CS = coordinate system
  real (*nucleationStressInFaultCS)[misc::numPaddedPoints][6];
  real dt = 0;
  real gNuc = 0;

  real (*a)[misc::numPaddedPoints];
  real (*sl0)[misc::numPaddedPoints];
  real (*stateVariable)[misc::numPaddedPoints];

  TPMethod tpMethod;

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

  RateAndStateBase(DRParameters& drParameters)
      : BaseFrictionLaw<RateAndStateBase<Derived, TPMethod>>::BaseFrictionLaw(drParameters),
        tpMethod(TPMethod(drParameters)) {}

  protected:
  public:
  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             TractionResults& tractionResults,
                             std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::numPaddedPoints>& strengthBuffer,
                             unsigned& ltsFace,
                             unsigned& timeIndex) {
    bool hasConverged = false;

    tpMethod.setInitialFluidPressure(ltsFace);
    // compute initial slip rate and reference values
    std::array<real, misc::numPaddedPoints> localSlipRate{0};
    std::array<real, misc::numPaddedPoints> stateVarReference{0};
    std::array<real, misc::numPaddedPoints> normalStress{0};
    std::array<real, misc::numPaddedPoints> absoluteShearStress{0};
    std::tie(absoluteShearStress, normalStress, stateVarReference, localSlipRate) =
        static_cast<Derived*>(this)->calcInitialVariables(
            faultStresses, stateVariableBuffer, tpMethod, timeIndex, ltsFace);
    // compute slip rates by solving non-linear system of equations (with newton)
    this->updateStateVariableIterative(hasConverged,
                                       stateVarReference,
                                       localSlipRate,
                                       stateVariableBuffer,
                                       tpMethod,
                                       normalStress,
                                       absoluteShearStress,
                                       faultStresses,
                                       timeIndex,
                                       ltsFace);

    // check for convergence
    if (!hasConverged) {
      static_cast<Derived*>(this)->executeIfNotConverged(stateVariableBuffer, ltsFace);
    }
    // compute final thermal pressure and normalStress
    tpMethod.calcFluidPressure(normalStress,
                               this->mu,
                               localSlipRate,
                               this->deltaT[timeIndex],
                               true,
                               timeIndex,
                               ltsFace);
    updateNormalStress(normalStress, faultStresses, tpMethod, timeIndex, ltsFace);
    // compute final slip rates and traction from median value of the iterative solution and initial
    // guess
    this->calcSlipRateAndTraction(stateVarReference,
                                  localSlipRate,
                                  stateVariableBuffer,
                                  normalStress,
                                  absoluteShearStress,
                                  faultStresses,
                                  tractionResults,
                                  timeIndex,
                                  ltsFace);
  }

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {
    // compute time increments (Gnuc)
    this->preCalcTime();

    // compute initial stress
    static_cast<Derived*>(this)->adjustInitialStress(ltsFace);
    // copy state variable from last time step
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      stateVariableBuffer[pointIndex] = this->stateVariable[ltsFace][pointIndex];
    }
  }

  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned ltsFace) {
    // resample state variables
    std::array<real, misc::numPaddedPoints> resampledStateVar =
        static_cast<Derived*>(this)->resampleStateVar(stateVariableBuffer, ltsFace);
    // write back state Variable to lts tree
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] = resampledStateVar[pointIndex];
    }
  }

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    auto* concreteLts = dynamic_cast<seissol::initializers::LTS_RateAndState*>(dynRup);
    a = layerData.var(concreteLts->rsA);
    nucleationStressInFaultCS = layerData.var(concreteLts->nucleationStressInFaultCS);
    sl0 = layerData.var(concreteLts->rsSl0);
    stateVariable = layerData.var(concreteLts->stateVariable);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    tpMethod.copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }
  /*
   * compute time increments (Gnuc)
   */
  void preCalcTime() {
    dt = 0;
    for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
      dt += this->deltaT[timeIndex];
    }
    gNuc = this->calcSmoothStepIncrement(this->mFullUpdateTime, dt);
  }

  /*
   * compute initial stress from nucleation stress
   */
  void adjustInitialStress(unsigned int ltsFace) {
    if (this->mFullUpdateTime <= this->drParameters.t0) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        for (unsigned i = 0; i < 6; i++) {
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
  std::tuple<std::array<real, misc::numPaddedPoints>,
             std::array<real, misc::numPaddedPoints>,
             std::array<real, misc::numPaddedPoints>,
             std::array<real, misc::numPaddedPoints>>
      calcInitialVariables(FaultStresses const& faultStresses,
                           std::array<real, misc::numPaddedPoints> const& localStateVariable,
                           TPMethod const& tpMethod,
                           unsigned int timeIndex,
                           unsigned int ltsFace) {
    // Careful, the state variable must always be corrected using stateVarZero and not
    // localStateVariable!
    std::array<real, misc::numPaddedPoints> stateVarZero;
    std::copy(localStateVariable.begin(), localStateVariable.end(), stateVarZero.begin());

    std::array<real, misc::numPaddedPoints> absoluteShearStress;
    std::array<real, misc::numPaddedPoints> normalStress;
    std::array<real, misc::numPaddedPoints> temporarySlipRate;

    updateNormalStress(normalStress, faultStresses, tpMethod, timeIndex, ltsFace);
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.xyStress[timeIndex][pointIndex];
      real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.xzStress[timeIndex][pointIndex];
      absoluteShearStress[pointIndex] = misc::magnitude(totalStressXY, totalStressXZ);


      // The following process is adapted from that described by Kaneko et al. (2008)
      this->slipRateMagnitude[ltsFace][pointIndex] = misc::magnitude(
          this->slipRate1[ltsFace][pointIndex], this->slipRate2[ltsFace][pointIndex]);
      this->slipRateMagnitude[ltsFace][pointIndex] =
          std::max(almostZero, this->slipRateMagnitude[ltsFace][pointIndex]);
      temporarySlipRate[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
    } // End of pointIndex-loop
    return {absoluteShearStress, normalStress, stateVarZero, temporarySlipRate};
  }

  void updateStateVariableIterative(
      bool& hasConverged,
      std::array<real, misc::numPaddedPoints> const& stateVarReference,
      std::array<real, misc::numPaddedPoints>& localSlipRate,
      std::array<real, misc::numPaddedPoints>& localStateVariable,
      TPMethod& tpMethod,
      std::array<real, misc::numPaddedPoints>& normalStress,
      std::array<real, misc::numPaddedPoints> const& absoluteShearStress,
      FaultStresses const& faultStresses,
      unsigned int timeIndex,
      unsigned int ltsFace) {
    std::array<real, misc::numPaddedPoints> testSlipRate{0};
    for (unsigned j = 0; j < numberStateVariableUpdates; j++) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        // fault strength using friction coefficient and fluid pressure from previous
        // timestep/iteration update state variable using sliprate from the previous time step
        localStateVariable[pointIndex] =
            static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                             ltsFace,
                                                             stateVarReference[pointIndex],
                                                             this->deltaT[timeIndex],
                                                             localSlipRate[pointIndex]);
      }
      tpMethod.calcFluidPressure(normalStress,
                                 this->mu,
                                 localSlipRate,
                                 this->deltaT[timeIndex],
                                 false,
                                 timeIndex,
                                 ltsFace);

      updateNormalStress(normalStress, faultStresses, tpMethod, timeIndex, ltsFace);

      // solve for new slip rate, applying the Newton-Raphson algorithm
      // effective normal stress including initial stresses and pore fluid pressure
      hasConverged = this->IterativelyInvertSR(
          ltsFace, localStateVariable, normalStress, absoluteShearStress, testSlipRate);

      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        // update local slip rate, now using V=(Vnew+Vold)/2
        // For the next SV update, use the mean slip rate between the initial guess and the one
        // found (Kaneko 2008, step 6)
        localSlipRate[pointIndex] =
            0.5 * (this->slipRateMagnitude[ltsFace][pointIndex] + fabs(testSlipRate[pointIndex]));

        // 4. solve again for Vnew
        this->slipRateMagnitude[ltsFace][pointIndex] = fabs(testSlipRate[pointIndex]);

        // update friction coefficient based on new state variable and slip rate
        this->mu[ltsFace][pointIndex] =
            static_cast<Derived*>(this)->updateMu(ltsFace,
                                                  pointIndex,
                                                  this->slipRateMagnitude[ltsFace][pointIndex],
                                                  localStateVariable[pointIndex]);
      } // End of pointIndex-loop
    }
  }

  void calcSlipRateAndTraction(std::array<real, misc::numPaddedPoints> const& stateVarReference,
                               std::array<real, misc::numPaddedPoints> const& localSlipRate,
                               std::array<real, misc::numPaddedPoints>& localStateVariable,
                               std::array<real, misc::numPaddedPoints> const& normalStress,
                               std::array<real, misc::numPaddedPoints> const& absoluteShearStress,
                               FaultStresses& faultStresses,
                               TractionResults& tractionResults,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      // SV from mean slip rate in tmp
      localStateVariable[pointIndex] =
          static_cast<Derived*>(this)->updateStateVariable(pointIndex,
                                                           ltsFace,
                                                           stateVarReference[pointIndex],
                                                           this->deltaT[timeIndex],
                                                           localSlipRate[pointIndex]);

      // update LocMu for next strength determination, only needed for last update
      this->mu[ltsFace][pointIndex] =
          static_cast<Derived*>(this)->updateMu(ltsFace,
                                                pointIndex,
                                                this->slipRateMagnitude[ltsFace][pointIndex],
                                                localStateVariable[pointIndex]);
      real strength = -this->mu[ltsFace][pointIndex] * normalStress[pointIndex];
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = this->initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.xyStress[timeIndex][pointIndex];
      real totalStressXZ = this->initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.xzStress[timeIndex][pointIndex];
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
      this->accumulatedSlipMagnitude[ltsFace][pointIndex] +=
          this->slipRateMagnitude[ltsFace][pointIndex] * this->deltaT[timeIndex];

      // Update slip rate (notice that locSlipRate(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate
      // caused by a free surface!)
      this->slipRate1[ltsFace][pointIndex] =
          -this->impAndEta[ltsFace].invEtaS *
          (this->tractionXY[ltsFace][pointIndex] - faultStresses.xyStress[timeIndex][pointIndex]);
      this->slipRate2[ltsFace][pointIndex] =
          -this->impAndEta[ltsFace].invEtaS *
          (this->tractionXZ[ltsFace][pointIndex] - faultStresses.xzStress[timeIndex][pointIndex]);

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
      tractionResults.xyTraction[timeIndex][pointIndex] = this->tractionXY[ltsFace][pointIndex];
      tractionResults.xzTraction[timeIndex][pointIndex] = this->tractionXZ[ltsFace][pointIndex];

      // update directional slip
      this->slip1[ltsFace][pointIndex] +=
          this->slipRate1[ltsFace][pointIndex] * this->deltaT[timeIndex];
      this->slip2[ltsFace][pointIndex] +=
          this->slipRate2[ltsFace][pointIndex] * this->deltaT[timeIndex];
    } // End of BndGP-loop
  }

  void saveDynamicStressOutput(unsigned int face) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

      if (this->ruptureTime[face][pointIndex] > 0.0 &&
          this->ruptureTime[face][pointIndex] <= this->mFullUpdateTime &&
          this->dynStressTimePending[pointIndex] &&
          this->mu[face][pointIndex] <=
              (this->drParameters.muW +
               0.05 * (this->drParameters.rsF0 - this->drParameters.muW))) {
        this->dynStressTime[face][pointIndex] = this->mFullUpdateTime;
        this->dynStressTimePending[face][pointIndex] = false;
      }
    }
  }

  bool IterativelyInvertSR(unsigned int ltsFace,
                           std::array<real, misc::numPaddedPoints> const& localStateVariable,
                           std::array<real, misc::numPaddedPoints> const& normalStress,
                           std::array<real, misc::numPaddedPoints> const& absoluteShearStress,
                           std::array<real, misc::numPaddedPoints>& slipRateTest) {

    // Note that we need double precision here, since single precision led to NaNs.
    double muF[misc::numPaddedPoints], dMuF[misc::numPaddedPoints];
    double nr[misc::numPaddedPoints], dNr[misc::numPaddedPoints];
    // double AlmostZero = 1e-45;
    bool hasConverged = false;

    // solve for Vnew = SR , applying the Newton-Raphson algorithm
    // SR fulfills g(SR)=f(SR)
    // -> find root of NR=f-g using a Newton-Raphson algorithm with dNr = d(NR)/d(SR)
    // SR_{i+1}=SR_i-( NR_i / dNR_i )
    //
    //        equalize:
    //         g = SR*MU/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
    //         f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
    //  where mu = friction coefficient, dependening on the RSF law used

    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      // first guess = SR value of the previous step
      slipRateTest[pointIndex] = this->slipRateMagnitude[ltsFace][pointIndex];
    }

    for (unsigned i = 0; i < numberSlipRateUpdates; i++) {
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

        // f = ( tmp2 * ABS(LocP+P_0)- ABS(S_0))*(S_0)/ABS(S_0)
        // g = slipRateTest * 1.0/(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) + ABS(ShTest)
        // for compiling reasons ASINH(X)=LOG(X+SQRT(X^2+1))

        // calculate friction coefficient
        muF[pointIndex] = static_cast<Derived*>(this)->updateMu(
            ltsFace, pointIndex, slipRateTest[pointIndex], localStateVariable[pointIndex]);
        dMuF[pointIndex] = static_cast<Derived*>(this)->updateMuDerivative(
            ltsFace, pointIndex, slipRateTest[pointIndex], localStateVariable[pointIndex]);
        nr[pointIndex] =
            -this->impAndEta[ltsFace].invEtaS * (fabs(normalStress[pointIndex]) * muF[pointIndex] -
                                                 absoluteShearStress[pointIndex]) -
            slipRateTest[pointIndex];
      }

      hasConverged = true;

      // max element of NR must be smaller than newtonTolerance
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
        if (fabs(nr[pointIndex]) >= newtonTolerance) {
          hasConverged = false;
          break;
        }
      }
      if (hasConverged) {
        return hasConverged;
      }
      for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {

        // derivative of NR
        dNr[pointIndex] = -this->impAndEta[ltsFace].invEtaS *
                              (fabs(normalStress[pointIndex]) * dMuF[pointIndex]) -
                          1.0;
        // ratio
        real tmp3 = nr[pointIndex] / dNr[pointIndex];

        // update slipRateTest
        slipRateTest[pointIndex] = std::max(almostZero, slipRateTest[pointIndex] - tmp3);
      }
    }
    return hasConverged;
  }

  void updateNormalStress(std::array<double, misc::numPaddedPoints>& normalStress, FaultStresses const& faultStresses, TPMethod tpMethod, size_t timeIndex, size_t ltsFace) {
    for (size_t pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      normalStress[pointIndex] = std::min(static_cast<real>(0.0),
                                          faultStresses.normalStress[timeIndex][pointIndex] +
                                          this->initialStressInFaultCS[ltsFace][pointIndex][0] -
                                          tpMethod.fluidPressure(ltsFace, pointIndex));
    }
  }
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_RATEANDSTATE_H
#include "LinearSlipWeakening.h"

void seissol::dr::friction_law::LinearSlipWeakeningLawFL2::calcStrengthHook(std::array<real, numOfPointsPadded> &Strength,
                                                                            FaultStresses &faultStresses,
                                                                            unsigned int iTimeGP,
                                                                            unsigned int ltsFace) {
  for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
    //-------------------------------------
    //calculate Fault Strength
    //fault strength (Uphoff eq 2.44) with addition cohesion term
    Strength[iBndGP] = cohesion[ltsFace][iBndGP] - mu[ltsFace][iBndGP] * std::min(initialStressInFaultCS[ltsFace][iBndGP][0] + faultStresses.NormalStressGP[iTimeGP][iBndGP], (real)0.0);
  }
}

void seissol::dr::friction_law::LinearSlipWeakeningLawFL2::calcStateVariableHook(std::array<real, numOfPointsPadded> &stateVariablePsi,
                                                                                 std::array<real, numOfPointsPadded> &outputSlip,
                                                                                 dynamicRupture::kernel::resampleParameter &resampleKrnl,
                                                                                 unsigned int iTimeGP,
                                                                                 unsigned int ltsFace) {
  real resampledSlipRate[numOfPointsPadded];
  resampleKrnl.resamplePar = SlipRateMagnitude[ltsFace];
  resampleKrnl.resampledPar = resampledSlipRate;  //output from execute

  //Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees of freedom
  //resampleMatrix first projects LocSR on the two-dimensional basis on the reference triangle with
  //degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial at the quadrature points
  resampleKrnl.execute();

  for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
    //-------------------------------------
    //integrate Sliprate To Get Slip = State Variable
    slip[ltsFace][iBndGP] = slip[ltsFace][iBndGP] + resampledSlipRate[iBndGP] * deltaT[iTimeGP];
    outputSlip[iBndGP] = outputSlip[iBndGP] + SlipRateMagnitude[ltsFace][iBndGP] * deltaT[iTimeGP];

    //-------------------------------------
    //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
    //actually slip is already the stateVariable for this FL, but to simplify the next equations we divide it here by d_C
    stateVariablePsi[iBndGP] = std::min(std::fabs(slip[ltsFace][iBndGP]) / d_c[ltsFace][iBndGP], (real)1.0);
  }
}

void seissol::dr::friction_law::LinearSlipWeakeningLawFL16::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                                                               seissol::initializers::DynamicRupture *dynRup,
                                                                               real fullUpdateTime) {
  //first copy all Variables from the Base Lts dynRup tree
  LinearSlipWeakeningLawFL2::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  //maybe change later to const_cast?
  seissol::initializers::LTS_LinearSlipWeakeningFL16 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningFL16 *>(dynRup);
  forced_rupture_time                           = layerData.var(ConcreteLts->forced_rupture_time);
  tn                                            = layerData.var(ConcreteLts->tn);
}

void seissol::dr::friction_law::LinearSlipWeakeningLawFL16:: setTimeHook(unsigned int ltsFace) {
  tn[ltsFace] = m_fullUpdateTime;
}

void seissol::dr::friction_law::LinearSlipWeakeningLawFL16::calcStateVariableHook(std::array<real, numOfPointsPadded> &stateVariablePsi,
                                                                                  std::array<real, numOfPointsPadded> &outputSlip,
                                                                                  dynamicRupture::kernel::resampleParameter &resampleKrnl,
                                                                                  unsigned int iTimeGP,
                                                                                  unsigned int ltsFace) {
  LinearSlipWeakeningLawFL2::calcStateVariableHook(stateVariablePsi, outputSlip, resampleKrnl, iTimeGP, ltsFace);
  tn[ltsFace] += deltaT[iTimeGP];

  for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
    real f2 = 0.0;
    if (m_Params->t_0 == 0) {
      if (tn[ltsFace] >= forced_rupture_time[ltsFace][iBndGP] ) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max((real)0.0, std::min( (real)1.0, (m_fullUpdateTime - forced_rupture_time[ltsFace][iBndGP] ) / m_Params->t_0));
    }
    stateVariablePsi[iBndGP] = std::max(stateVariablePsi[iBndGP], f2);
  }
}

void seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6::calcStrengthHook(std::array<real, numOfPointsPadded> &Strength,
                                                                                     FaultStresses &faultStresses,
                                                                                     unsigned int iTimeGP,
                                                                                     unsigned int ltsFace) {
  std::array<real, numOfPointsPadded> LocSlipRate;
  std::array<real, numOfPointsPadded> sigma;

  for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
    //  modify strength according to prakash clifton
    // literature e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
    LocSlipRate[iBndGP] = std::sqrt(slipRateStrike[ltsFace][iBndGP] * slipRateStrike[ltsFace][iBndGP] + slipRateDip[ltsFace][iBndGP] * slipRateDip[ltsFace][iBndGP]);
    sigma[iBndGP] = faultStresses.NormalStressGP[iTimeGP][iBndGP] + initialStressInFaultCS[ltsFace][iBndGP][0];
    prak_clif_mod(strengthData[ltsFace][iBndGP], sigma[iBndGP], LocSlipRate[iBndGP], mu[ltsFace][iBndGP], deltaT[iTimeGP]);

    //TODO: add this line to make the FL6 actually functional: (this line is also missing in the master branch)
    //Strength[iBndGP] = strengthData[ltsFace][iBndGP];
  }
}

void seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6::calcStateVariableHook(std::array<real, numOfPointsPadded> &stateVariablePsi,
                                                                                           std::array<real, numOfPointsPadded> &outputSlip,
                                                                                           dynamicRupture::kernel::resampleParameter &resampleKrnl,
                                                                                           unsigned int iTimeGP,
                                                                                           unsigned int ltsFace) {
  for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
    slip[ltsFace][iBndGP] = slip[ltsFace][iBndGP] + SlipRateMagnitude[ltsFace][iBndGP] * deltaT[iTimeGP];
    outputSlip[iBndGP] = slip[ltsFace][iBndGP];

    //-------------------------------------
    //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
    //actually slip is already the stateVariable for this FL, but to simplify the next equations we divide it here by d_C
    stateVariablePsi[iBndGP] = std::min(std::fabs(slip[ltsFace][iBndGP]) / d_c[ltsFace][iBndGP], (real)1.0);
  }
}


void seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6::copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                                                                                        seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) {
//first copy all Variables from the Base Lts dynRup tree
LinearSlipWeakeningLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
//maybe change later to const_cast?
seissol::initializers::LTS_LinearBimaterialFL6 *ConcreteLts = dynamic_cast<seissol::initializers::LTS_LinearBimaterialFL6 *>(dynRup);
strengthData = layerData.var(ConcreteLts->strengthData);
}

/*
 * calculates strength
 */
void seissol::dr::friction_law::LinearSlipWeakeningLawBimaterialFL6::prak_clif_mod(real &strength, real &sigma, real &LocSlipRate, real &mu, real &dt){
  real expterm;
  expterm = std::exp(-(std::abs(LocSlipRate) + m_Params->v_star)*dt/ m_Params->prakash_length);
  strength =  strength*expterm - std::max((real)0.0,-mu*sigma)*(expterm-1.0);
}
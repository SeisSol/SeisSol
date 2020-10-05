//
// Created by adrian on 17.09.20.
//

#ifndef SEISSOL_DR_SOLVER_LINEAR_H
#define SEISSOL_DR_SOLVER_LINEAR_H


#include "DR_solver_base.h"

#include "DR_math.h"
#include <yaml-cpp/yaml.h>

namespace seissol {
  namespace dr {
    namespace fr_law {
      class LinearSlipWeakeningSolverFL2;  //linear slip weakening
      class Solver_FL_16; //Linear slip weakening forced time rapture
      class LinearSlipWeakSolverBimaterialFL6;
    }
  }
}


/*
*     !> friction case 16,17
*     !> Specific conditions for SCEC TPV16/17
*     !> basically, introduction of a time dependent forced rupture
*/
class seissol::dr::fr_law::LinearSlipWeakeningSolverFL2 : public seissol::dr::fr_law::BaseFrictionSolver {
protected:
  //parameter for insta_healing
  //TODO: make this parameter better accessible?
  real u_0  = 10e-14; //slip rate is considered as being zero for instaneous healing
  real                    (*d_c)[numOfPointsPadded];
  real                    (*mu_S)[numOfPointsPadded];
  real                    (*mu_D)[numOfPointsPadded];
  bool                    (*DS)[numOfPointsPadded];
  real                    (*dynStress_time)[numOfPointsPadded];
  //only for FL16:
  real                    (*forced_rupture_time)[numOfPointsPadded];
  real                    *tn;


  //Hook for Factory_FL_16
  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real fullUpdateTime, unsigned int iBndGP, unsigned int ltsFace  ) {
    //!Do nothing
  }

  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup) override {
    //first copy all Variables from the Base Lts dynRup tree
    BaseFrictionSolver::copyLtsTreeToLocal(layerData, dynRup);
    //TODO: change later to const_cast
    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);
    d_c                                           = layerData.var(ConcreteLts->d_c);
    mu_S                                          = layerData.var(ConcreteLts->mu_S);
    mu_D                                          = layerData.var(ConcreteLts->mu_D);
    DS                                            = layerData.var(ConcreteLts->DS);
    averaged_Slip                                 = layerData.var(ConcreteLts->averaged_Slip);
    dynStress_time                                = layerData.var(ConcreteLts->dynStress_time);

    //only for FL16:
    forced_rupture_time                           = layerData.var(ConcreteLts->forced_rupture_time);
    tn                                            = layerData.var(ConcreteLts->tn);
  }

  /*
   *
   */
  void calcSlipRateAndTraction(
      FaultStresses &faultStresses,
      real LocSlipRate[seissol::tensor::resamplePar::size()],
      real DeltaT,
      unsigned int iTimeGP,
      unsigned int face
  ){
    std::array<real, numOfPointsPadded> TotalShearStressYZ;
    std::array<real, numOfPointsPadded> Strength;

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      //-------------------------------------
      //calculate Fault Strength
      //fault strength (Uphoff eq 2.44)
      Strength[iBndGP] = cohesion[face][iBndGP] - mu[face][iBndGP] * std::min(initialStressInFaultCS[face][iBndGP][0] + faultStresses.NorStressGP[iTimeGP][iBndGP], 0.0);

      //-------------------------------------
      //calculate TotalShearStress in Y and Z direction
      TotalShearStressYZ[iBndGP] = std::sqrt(
          seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP], 2) +
          seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP], 2));

      //-------------------------------------
      // calculate SlipRates
      LocSlipRate[iBndGP] = std::max(0.0, (TotalShearStressYZ[iBndGP] - Strength[iBndGP]) / impAndEta[face].eta_s);
      slipRate1[face][iBndGP] = LocSlipRate[iBndGP] * (initialStressInFaultCS[face][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP]) /
                                (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));
      slipRate2[face][iBndGP]  = LocSlipRate[iBndGP] * (initialStressInFaultCS[face][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP]) /
                                 (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));

      //-------------------------------------
      //calculateTraction
      faultStresses.TractionGP_XY[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP] - impAndEta[face].eta_s * slipRate1[face][iBndGP];
      faultStresses.TractionGP_XZ[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP] - impAndEta[face].eta_s * slipRate2[face][iBndGP];
      tracXY[face][iBndGP] = faultStresses.TractionGP_XY[iTimeGP][iBndGP];
      tracXZ[face][iBndGP] = faultStresses.TractionGP_XY[iTimeGP][iBndGP];

      //-------------------------------------
      //update Directional Slip
      slip1[face][iBndGP] = slip1[face][iBndGP] + slipRate1[face][iBndGP] * DeltaT;
      slip2[face][iBndGP] = slip2[face][iBndGP] + slipRate2[face][iBndGP] * DeltaT;
    }
  }

/*
 * tmpSlip output
 */
  void calcStateVariableAndFrictionFunc(
      std::array<real, numOfPointsPadded> &tmpSlip,
      real LocSlipRate[seissol::tensor::resamplePar::size()],
      dynamicRupture::kernel::resampleParameter &resampleKrnl,
      real fullUpdateTime,
      real DeltaT,
      unsigned int ltsFace
  ){
    std::array<real, numOfPointsPadded> stateVariablePsi;
    real resampledSlipRate[seissol::tensor::resamplePar::size()];
    resampleKrnl.resamplePar = LocSlipRate;
    resampleKrnl.resampledPar = resampledSlipRate;  //output from execute

    //Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees of freedom
    //resampleMatrix first projects LocSR on the two-dimensional basis on the reference triangle with
    //degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial at the quadrature points
    resampleKrnl.execute();



    //TODO: does not work with padded Points bc of resampleMatrix is not padded
    for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
      //-------------------------------------
      //integrate Sliprate To Get Slip = State Variable
      slip[ltsFace][iBndGP] = slip[ltsFace][iBndGP] + resampledSlipRate[iBndGP] * DeltaT;
      tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSlipRate[iBndGP] * DeltaT;

      //-------------------------------------
      //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
      //actually slip is already the stateVariable for this FL, but to simplify the next equations we divide it here by d_C
      stateVariablePsi[iBndGP] = std::min(std::fabs(slip[ltsFace][iBndGP]) / d_c[ltsFace][iBndGP], 1.0);

      //-------------------------------------
      //hook for Factory_FL_16: may calculate a different value for the state variable Psi
      hookCalcStateVariable(stateVariablePsi, fullUpdateTime, iBndGP, ltsFace); //for Factory_FL_16

      //-------------------------------------
      //Carsten Thesis: Eq. 2.45
      //evaluate friction law: updated mu -> friction law
      mu[ltsFace][iBndGP] = mu_S[ltsFace][iBndGP] - (mu_S[ltsFace][iBndGP] - mu_D[ltsFace][iBndGP]) * stateVariablePsi[iBndGP];

      //-------------------------------------
      //instantaneous healing option Reset Mu and Slip
      if (m_Params->IsInstaHealingOn == true) {
        if (LocSlipRate[iBndGP] < u_0) {
          mu[ltsFace][iBndGP] = mu_S[ltsFace][iBndGP];
          slip[ltsFace][iBndGP] = 0.0;
        }
      }
    }//end of iBndGP-loop
  }//end of function calcStateVariableAndFrictionFunc


  //output time when shear stress is equal to the dynamic stress after rupture arrived
  //currently only for linear slip weakening
  void outputDynamicStress(
      real fullUpdateTime,
      unsigned int face
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      if (rupture_time[face][iBndGP] > 0.0 &&
          rupture_time[face][iBndGP] <= fullUpdateTime &&
          DS[iBndGP] &&
          std::fabs(slip[face][iBndGP]) >= d_c[face][iBndGP]) {
        dynStress_time[face][iBndGP] = fullUpdateTime;
        DS[face][iBndGP] = false;
      }
    }
  }

public:
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER],
                        real DeltaT[CONVERGENCE_ORDER]) override {

    copyLtsTreeToLocal(layerData, dynRup);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      //initialize struct for in/outputs stresses
      FaultStresses faultStresses{};

      //declare local variables
      std::array<real, numOfPointsPadded> tmpSlip{0};
      real LocSlipRate[seissol::tensor::resamplePar::size()];
      dynamicRupture::kernel::resampleParameter resampleKrnl;
      resampleKrnl.resampleM = init::resample::Values;
      //tn only needed for FL=16
      tn[ltsFace] = fullUpdateTime;

      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        tn[ltsFace] += DeltaT[iTimeGP];
        calcSlipRateAndTraction(faultStresses, LocSlipRate, DeltaT[iTimeGP], iTimeGP , ltsFace);

        calcStateVariableAndFrictionFunc(tmpSlip, LocSlipRate, resampleKrnl, fullUpdateTime, DeltaT[iTimeGP], ltsFace);
      }//End of iTimeGP-Loop

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(LocSlipRate, fullUpdateTime, ltsFace);

      //output time when shear stress is equal to the dynamic stress after rupture arrived
      //currently only for linear slip weakening
      outputDynamicStress(fullUpdateTime, ltsFace);

      //output peak slip rate
      calcPeakSlipRate(LocSlipRate, ltsFace);

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      calcAverageSlip(tmpSlip, ltsFace);

      postcomputeImposedStateFromNewStress(
          QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace],
          faultStresses, timeWeights, ltsFace);
    }//End of Loop over Faces
  }//End of Function evaluate

};//End of Class



class seissol::dr::fr_law::Solver_FL_16 : public seissol::dr::fr_law::LinearSlipWeakeningSolverFL2 {
public:

  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real fullUpdateTime,  unsigned int iBndGP,  unsigned int ltsFace) override{
    real f2 = 0.0;

    if (m_Params->t_0 == 0) {
      if (tn[ltsFace] >= forced_rupture_time[ltsFace][iBndGP] ) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max(0.0, std::min( 1.0, (fullUpdateTime - forced_rupture_time[ltsFace][iBndGP] ) / m_Params->t_0));
    }
    stateVariablePsi[iBndGP] = std::max(stateVariablePsi[iBndGP], f2);
  }
};


class seissol::dr::fr_law::LinearSlipWeakSolverBimaterialFL6 : public seissol::dr::fr_law::LinearSlipWeakeningSolverFL2 {
protected:
  //Attributes
  real  (*strengthData)[numOfPointsPadded];

  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup) override {
    //first copy all Variables from the Base Lts dynRup tree
    LinearSlipWeakeningSolverFL2::copyLtsTreeToLocal(layerData, dynRup);
    //TODO: change later to const_cast
    seissol::initializers::DR_FL_6 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_6 *>(dynRup);
    strengthData = layerData.var(ConcreteLts->strengthData);
  }

  /*
   * calculates strength
   */
  void prak_clif_mod(real &strength, real &sigma, real &LocSlipRate, real &mu, real &dt){
    real expterm;
    expterm = std::exp(-(std::abs(LocSlipRate) + m_Params->v_star)*dt/ m_Params->prakash_length);
    strength =  strength* expterm - std::max(0.0,-mu*sigma)*(expterm-1.0);
  }

public:
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER],
                        real DeltaT[CONVERGENCE_ORDER]) override {
    copyLtsTreeToLocal(layerData, dynRup);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      //initialize struct for in/outputs stresses
      FaultStresses faultStresses{};

      //declare local variables
      real LocSlipRate[seissol::tensor::resamplePar::size()];

      //compute stresses from Qinterpolated
      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      real sigma;
      real ShTest;
      real Strength;

      for(int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++){ //loop over time steps
        //TODO: test if works padded??
        for(int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++){ //loop over all points

          //  modify strength according to prakash clifton
          LocSlipRate[iBndGP] = std::sqrt(slipRate1[ltsFace][iBndGP]*slipRate1[ltsFace][iBndGP] + slipRate2[ltsFace][iBndGP]*slipRate2[ltsFace][iBndGP]);
          sigma = faultStresses.NorStressGP[iTimeGP][iBndGP]+initialStressInFaultCS[ltsFace][iBndGP][0];
          prak_clif_mod(strengthData[ltsFace][iBndGP], sigma, LocSlipRate[iBndGP], mu[ltsFace][iBndGP], DeltaT[iTimeGP]);
          ShTest = sqrt(seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP], 2) + seissol::dr::aux::power(initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP],2) );

          if(ShTest > Strength){
            // 1 evaluate friction
            faultStresses.TractionGP_XY[iTimeGP][iBndGP] = ((initialStressInFaultCS[ltsFace][iBndGP][3] + faultStresses.XYStressGP[iTimeGP][iBndGP])/ShTest)*Strength;
            faultStresses.TractionGP_XZ[iTimeGP][iBndGP] = ((initialStressInFaultCS[ltsFace][iBndGP][5] + faultStresses.XZStressGP[iTimeGP][iBndGP])/ShTest)*Strength;

            // 2 update stress change
            faultStresses.TractionGP_XY[iTimeGP][iBndGP] -= initialStressInFaultCS[ltsFace][iBndGP][3];
            faultStresses.TractionGP_XZ[iTimeGP][iBndGP] -= initialStressInFaultCS[ltsFace][iBndGP][5];
          }else{
            faultStresses.TractionGP_XY[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP];
            faultStresses.TractionGP_XZ[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP];
          }

          //!Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
          slipRate1[ltsFace][iBndGP]     = -impAndEta[ltsFace].inv_eta_s*(faultStresses.TractionGP_XY[iTimeGP][iBndGP]-faultStresses.XYStressGP[iTimeGP][iBndGP]);
          slipRate2[ltsFace][iBndGP]     = -impAndEta[ltsFace].inv_eta_s*(faultStresses.TractionGP_XZ[iTimeGP][iBndGP]-faultStresses.XZStressGP[iTimeGP][iBndGP]);
          LocSlipRate[iBndGP]      = sqrt(slipRate1[ltsFace][iBndGP]*slipRate1[ltsFace][iBndGP] + slipRate2[ltsFace][iBndGP]*slipRate2[ltsFace][iBndGP]);
          // Update slip
          slip1[ltsFace][iBndGP] = slip1[ltsFace][iBndGP] + slipRate1[ltsFace][iBndGP]*DeltaT[iTimeGP];
          slip2[ltsFace][iBndGP] = slip2[ltsFace][iBndGP] + slipRate2[ltsFace][iBndGP]*DeltaT[iTimeGP];
          slip[ltsFace][iBndGP] = slip[ltsFace][iBndGP] + LocSlipRate[iBndGP]*DeltaT[iTimeGP];

          //friction function
          if(abs(slip[ltsFace][iBndGP]) < d_c[ltsFace][iBndGP]){
            mu[ltsFace][iBndGP] = mu_S[ltsFace][iBndGP] - (mu_S[ltsFace][iBndGP]-mu_D[ltsFace][iBndGP])/d_c[ltsFace][iBndGP]*abs(slip[ltsFace][iBndGP]);
          }else{
            mu[ltsFace][iBndGP] = mu_D[ltsFace][iBndGP];
          }

          // instantaneous healing
          if( m_Params->IsInstaHealingOn == true){
            if(LocSlipRate[iBndGP] < u_0){
              mu[ltsFace][iBndGP] = mu_S[ltsFace][iBndGP];
              // reset slip history for LSW
              slip[ltsFace][iBndGP] = 0.0; //0.0D0
            }
          }
        }//End loop over all points
      } //End of time loop

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(LocSlipRate, fullUpdateTime, ltsFace);

      //output peak slip rate
      calcPeakSlipRate(LocSlipRate, ltsFace);

      //output time when shear stress is equal to the dynamic stress after rupture arrived
      //currently only for linear slip weakening
      outputDynamicStress(fullUpdateTime, ltsFace);

      //save stresses in imposedState
      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);

    }//End of Loop over Faces

  }//End of Function evaluate
};


#endif //SEISSOL_DR_SOLVER_LINEAR_H

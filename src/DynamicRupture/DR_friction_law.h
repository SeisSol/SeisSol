//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_FRICTION_LAW_H
#define SEISSOL_DR_FRICTION_LAW_H

#include <c++/8.3.0/iostream>
#include "DR_math.h"


namespace seissol {
  namespace dr {
    namespace fr_law {
      class Base;
      class FL_2;
      class FL_3; //aging law
      class FL_4; //slip law
      class FL_16;
      class FL_33;
    }
  }
}


class seissol::dr::fr_law::Base {

public:
  //TODO: rename e.g. BaseSolverFL
  virtual ~Base() {}

protected:
  static constexpr int numberOfPoints =  tensor::QInterpolated::Shape[0];// DISC%Galerkin%nBndGP
  //TODO: is init::QInterpolated::Start[0] always 0?
  //assert(init::QInterpolated::Start[0] == 0);
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
  ImpedancesAndEta*                     impAndEta;
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus;
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus;
  real                    (*initialStressInFaultCS)[numOfPointsPadded][6];
  real                    (*cohesion)[numOfPointsPadded];
  real                    (*mu)[numOfPointsPadded];
  real                    (*slip)[numOfPointsPadded];
  real                    (*slip1)[numOfPointsPadded];
  real                    (*slip2)[numOfPointsPadded];
  real                    (*slipRate1)[numOfPointsPadded];
  real                    (*slipRate2)[numOfPointsPadded];
  real                    (*rupture_time)[numOfPointsPadded];
  bool                    (*RF)[numOfPointsPadded];
  real                    (*peakSR)[numOfPointsPadded];
  //TODO: merge TractionGP_XY and tracXY in one variable
  real                    (*tracXY)[numOfPointsPadded];
  real                    (*tracXZ)[numOfPointsPadded];
  real                    (*imposedStatePlus)[tensor::QInterpolated::size()];
  real                    (*imposedStateMinus)[tensor::QInterpolated::size()];

  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  virtual void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup){
    impAndEta                                     = layerData.var(dynRup->impAndEta);
    waveSpeedsPlus                                = layerData.var(dynRup->waveSpeedsPlus);
    waveSpeedsMinus                               = layerData.var(dynRup->waveSpeedsMinus);
    initialStressInFaultCS                        = layerData.var(dynRup->initialStressInFaultCS);
    cohesion                                      = layerData.var(dynRup->cohesion);
    mu                                            = layerData.var(dynRup->mu);
    slip                                          = layerData.var(dynRup->slip);
    slip1                                         = layerData.var(dynRup->slip1);
    slip2                                         = layerData.var(dynRup->slip2);
    slipRate1                                     = layerData.var(dynRup->slipRate1);
    slipRate2                                     = layerData.var(dynRup->slipRate2);
    rupture_time                                  = layerData.var(dynRup->rupture_time);
    RF                                            = layerData.var(dynRup->RF);
    peakSR                                        = layerData.var(dynRup->peakSR);
    tracXY                                        = layerData.var(dynRup->tracXY);
    tracXZ                                        = layerData.var(dynRup->tracXZ);
    imposedStatePlus                              = layerData.var(dynRup->imposedStatePlus);
    imposedStateMinus                             = layerData.var(dynRup->imposedStateMinus);
  }

  /*
   * output:
   * NorStressGP, XYStressGP, XZStressGP
   *
   * input:
   * QInterpolatedPlus, QInterpolatedMinus, eta_p, Zp, Zp_neig, eta_s, Zs, Zs_neig
   *
   * Calculate stress from jump of plus and minus side
   */
  void precomputeStressFromQInterpolated(
    real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
    real XYStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
    real XZStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
    real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    unsigned int face
    ){

    for(int j = 0; j < CONVERGENCE_ORDER; j++){
      auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
      auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      //TODO: does QInterpolatedMinusView work with padded access?
      for(int i = 0; i < numberOfPoints; i++){
        //Carsten Uphoff Thesis: EQ.: 4.53
        NorStressGP[j][i] = impAndEta[face].eta_p * (QInterpolatedMinusView(i,6) - QInterpolatedPlusView(i,6) + QInterpolatedPlusView(i,0) / impAndEta[face].Zp + QInterpolatedMinusView(i,0) / impAndEta[face].Zp_neig);
        XYStressGP[j][i]  = impAndEta[face].eta_s * (QInterpolatedMinusView(i,7) - QInterpolatedPlusView(i,7) + QInterpolatedPlusView(i,3) / impAndEta[face].Zs + QInterpolatedMinusView(i,3) / impAndEta[face].Zs_neig);
        XZStressGP[j][i] = impAndEta[face].eta_s * (QInterpolatedMinusView(i,8) - QInterpolatedPlusView(i,8) + QInterpolatedPlusView(i,5) / impAndEta[face].Zs + QInterpolatedMinusView(i,5) / impAndEta[face].Zs_neig);
      }
    }
    //TODO: is this assert really needed?
    static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],"Different number of quadrature points?");
  }//End of precompute Function


  /*
   * Output: imposedStatePlus, imposedStateMinus
   *
   * Integrate over all Time points with the time weights and calculate the traction vor each side according to
   * Carsten Uphoff Thesis: EQ.: 4.60
   */
  void postcomputeImposedStateFromNewStress(
      real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
      real TractionGP_XY[CONVERGENCE_ORDER][numOfPointsPadded],
      real TractionGP_XZ[CONVERGENCE_ORDER][numOfPointsPadded],
      real timeWeights[CONVERGENCE_ORDER],
      unsigned int face
      ){
    auto imposedStatePlusView = init::QInterpolated::view::create(imposedStatePlus[face]);
    auto imposedStateMinusView = init::QInterpolated::view::create(imposedStateMinus[face]);
    //initialize to 0
    imposedStateMinusView.setZero();
    imposedStatePlusView.setZero();

    for (int j = 0; j < CONVERGENCE_ORDER; j++) {
      auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
      auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      for (int i = 0; i < numberOfPoints; i++) {
        //Carsten Uphoff Thesis: EQ.: 4.60
        imposedStateMinusView(i, 0) += timeWeights[j] * NorStressGP[j][i];
        imposedStateMinusView(i, 3) += timeWeights[j] * TractionGP_XY[j][i];
        imposedStateMinusView(i, 5) += timeWeights[j] * TractionGP_XZ[j][i];
        imposedStateMinusView(i, 6) += timeWeights[j] * (QInterpolatedMinusView(i, 6) -  (NorStressGP[j][i] - QInterpolatedMinusView(i, 0))/  impAndEta[face].Zp_neig);
        imposedStateMinusView(i, 7) += timeWeights[j] * (QInterpolatedMinusView(i, 7) -  (TractionGP_XY[j][i] - QInterpolatedMinusView(i, 3))/  impAndEta[face].Zs_neig);
        imposedStateMinusView(i, 8) += timeWeights[j] * (QInterpolatedMinusView(i, 8) -  (TractionGP_XZ[j][i] - QInterpolatedMinusView(i, 5))/  impAndEta[face].Zs_neig);

        imposedStatePlusView(i, 0) += timeWeights[j] * NorStressGP[j][i];
        imposedStatePlusView(i, 3) += timeWeights[j] * TractionGP_XY[j][i];
        imposedStatePlusView(i, 5) += timeWeights[j] * TractionGP_XZ[j][i];
        imposedStatePlusView(i, 6) += timeWeights[j] * (QInterpolatedPlusView(i, 6) +  (NorStressGP[j][i] - QInterpolatedPlusView(i, 0)) /  impAndEta[face].Zp);
        imposedStatePlusView(i, 7) += timeWeights[j] * (QInterpolatedPlusView(i, 7) +  (TractionGP_XY[j][i] - QInterpolatedPlusView(i, 3)) / impAndEta[face].Zs);
        imposedStatePlusView(i, 8) += timeWeights[j] * (QInterpolatedPlusView(i, 8) +  (TractionGP_XZ[j][i] - QInterpolatedPlusView(i, 5)) / impAndEta[face].Zs);
      } //End numberOfPoints-loop
    } //End CONVERGENCE_ORDER-loop

  }

  // output rupture front
  // outside of iTimeGP loop in order to safe an 'if' in a loop
  // this way, no subtimestep resolution possible
  void outputRuptureFront(
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real fullUpdateTime,
      unsigned int face
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (RF[face][iBndGP] && LocSlipRate[iBndGP] > 0.001) {
        rupture_time[face][iBndGP] = fullUpdateTime;
        RF[face][iBndGP] = false;
      }
    }
  }


  void calcPeakSlipRate(
      std::array<real, numOfPointsPadded> &LocSlipRate,
      unsigned int face){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (LocSlipRate[iBndGP] > peakSR[face][iBndGP]) {
        peakSR[face][iBndGP] = LocSlipRate[iBndGP];
      }
    }
  }

public:
  //TODO: change arguments to "const double& ref" for all arguments that are not changed (only input)
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                         seissol::initializers::DynamicRupture *dynRup,
                         real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                         real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                         real fullUpdateTime,
                         real timeWeights[CONVERGENCE_ORDER],
                         real DeltaT[CONVERGENCE_ORDER]) = 0;
};


/*
*     !> friction case 16,17
*     !> Specific conditions for SCEC TPV16/17
*     !> basically, introduction of a time dependent forced rupture
*/
class seissol::dr::fr_law::FL_2 : public seissol::dr::fr_law::Base {
protected:
  //parameter for insta_healing
  //TODO: make this parameter better accessible?
  real u_0  = 10e-14; //slip rate is considered as being zero for instaneous healing
  yateto::DenseTensorView<2,double,unsigned> resampleMatrixView = init::resample::view::create(const_cast<double *>(init::resample::Values));
  real                    (*d_c)[numOfPointsPadded];
  real                    (*mu_S)[numOfPointsPadded];
  real                    (*mu_D)[numOfPointsPadded];
  bool                    (*inst_healing);
  bool                    (*magnitude_out);
  bool                    (*DS)[numOfPointsPadded];
  real                    (*averaged_Slip);
  real                    (*dynStress_time)[numOfPointsPadded];
  //only for FL16:
  real                    (*forced_rupture_time)[numOfPointsPadded];
  real                    *lts_t_0;

protected:
  //Hook for FL_16
  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real tn, real fullUpdateTime, unsigned int iBndGP, unsigned int face  ) {
    //!Do nothing
  }

  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup) override {
    //first copy all Variables from the Base Lts dynRup tree
    Base::copyLtsTreeToLocal(layerData, dynRup);
    //TODO: change later to const_cast
    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);
    d_c                                           = layerData.var(ConcreteLts->d_c);
    mu_S                                          = layerData.var(ConcreteLts->mu_S);
    mu_D                                          = layerData.var(ConcreteLts->mu_D);
    inst_healing                                  = layerData.var(ConcreteLts->inst_healing);
    magnitude_out                                 = layerData.var(ConcreteLts->magnitude_out);

    DS                                            = layerData.var(ConcreteLts->DS);
    averaged_Slip                                 = layerData.var(ConcreteLts->averaged_Slip);

    dynStress_time                                = layerData.var(ConcreteLts->dynStress_time);

    //only for FL16:
    forced_rupture_time                           = layerData.var(ConcreteLts->forced_rupture_time);
    lts_t_0                                       = layerData.var(ConcreteLts->t_0);
  }

  /*
   *
   */
  void calcSlipRateAndTraction(
      real TractionGP_XY[numOfPointsPadded],
      real TractionGP_XZ[numOfPointsPadded],
      real NorStressGP[numOfPointsPadded],
      real XYStressGP[numOfPointsPadded],
      real XZStressGP[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real DeltaT,
      unsigned int face
      ){
    std::array<real, numOfPointsPadded> TotalShearStressYZ;
    std::array<real, numOfPointsPadded> Strength;

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
    //-------------------------------------
    //calculate Fault Strength
      //fault strength (Uphoff eq 2.44)
      Strength[iBndGP] = cohesion[face][iBndGP] - mu[face][iBndGP] * std::min(initialStressInFaultCS[face][iBndGP][0] + NorStressGP[iBndGP], 0.0);

    //-------------------------------------
    //calculate TotalShearStress in Y and Z direction
      TotalShearStressYZ[iBndGP] = std::sqrt(
          seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][3] + XYStressGP[iBndGP], 2) +
          seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][5] + XZStressGP[iBndGP], 2));

    //-------------------------------------
    // calculate SlipRates
      LocSlipRate[iBndGP] = std::max(0.0, (TotalShearStressYZ[iBndGP] - Strength[iBndGP]) / impAndEta[face].eta_s);
      slipRate1[face][iBndGP] = LocSlipRate[iBndGP] * (initialStressInFaultCS[face][iBndGP][3] + XYStressGP[iBndGP]) /
                       (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));
      slipRate2[face][iBndGP]  = LocSlipRate[iBndGP] * (initialStressInFaultCS[face][iBndGP][5] + XZStressGP[iBndGP]) /
                       (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));

    //-------------------------------------
    //calculateTraction
      TractionGP_XY[iBndGP] = XYStressGP[iBndGP] - impAndEta[face].eta_s * slipRate1[face][iBndGP];
      TractionGP_XZ[iBndGP] = XZStressGP[iBndGP] - impAndEta[face].eta_s * slipRate2[face][iBndGP];
      tracXY[face][iBndGP] = TractionGP_XY[iBndGP];
      tracXZ[face][iBndGP] = TractionGP_XY[iBndGP];

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
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real tn,
      real fullUpdateTime,
      real DeltaT,
      unsigned int face
      ){
    std::array<real, numOfPointsPadded> resampledSlipRate{0};
    std::array<real, numOfPointsPadded> stateVariablePsi;

    for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
      //-------------------------------------
      //integrate Sliprate To Get Slip = State Variable
      //TODO: does not work with padded Points bc of resampleMatrix is not padded
      for (int j = 0; j < numberOfPoints; j++) {
        //Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees of freedom
        //resampleMatrix first projects LocSR on the two-dimensional basis on the reference triangle with
        //degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial at the quadrature points
        resampledSlipRate[iBndGP] += resampleMatrixView(iBndGP, j) * LocSlipRate[j];
      }
      slip[face][iBndGP] = slip[face][iBndGP] + resampledSlipRate[iBndGP] * DeltaT;
      tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSlipRate[iBndGP] * DeltaT;

      //-------------------------------------
      //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
      //actually slip is already the stateVariable for this FL, but to simplify the next equations we divide it here by d_C
      stateVariablePsi[iBndGP] = std::min(std::abs(slip[face][iBndGP]) / d_c[face][iBndGP], 1.0);

      //-------------------------------------
      //hook for FL_16: may calculate a different value for the state variable Psi
      hookCalcStateVariable(stateVariablePsi, tn, fullUpdateTime, iBndGP, face); //for FL_16

      //-------------------------------------
      //Carsten Thesis: Eq. 2.45
      //evaluate friction law: updated mu -> friction law
      mu[face][iBndGP] = mu_S[face][iBndGP] - (mu_S[face][iBndGP] - mu_D[face][iBndGP]) * stateVariablePsi[iBndGP];

      //-------------------------------------
      //instantaneous healing option Reset Mu and Slip
      if (inst_healing[face] == true) {
        if (LocSlipRate[iBndGP] < u_0) {
          mu[face][iBndGP] = mu_S[face][iBndGP];
          slip[face][iBndGP] = 0.0;
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
          std::abs(slip[face][iBndGP]) >= d_c[face][iBndGP]) {
        dynStress_time[face][iBndGP] = fullUpdateTime;
        DS[face][iBndGP] = false;
      }
    }
  }

  //---compute and store slip to determine the magnitude of an earthquake ---
  //    to this end, here the slip is computed and averaged per element
  //    in calc_seissol.f90 this value will be multiplied by the element surface
  //    and an output happened once at the end of the simulation
  void calcAverageSlip(
      std::array<real, numOfPointsPadded> &tmpSlip,
      unsigned int face
  ){
    real sum_tmpSlip = 0;
    if (magnitude_out[face]) {
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++)
        sum_tmpSlip += tmpSlip[iBndGP];
      averaged_Slip[face] = averaged_Slip[face] + sum_tmpSlip / numberOfPoints;
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
    for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {

      //TODO: merge TractionGP_XY and tracXY in one variable
      real TractionGP_XY[CONVERGENCE_ORDER][numOfPointsPadded] = {{}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real TractionGP_XZ[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};// OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
      real XYStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
      real XZStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};

      //declare local variables
      std::array<real, numOfPointsPadded> tmpSlip{0};
      std::array<real, numOfPointsPadded> LocSlipRate;
      //tn only needed for FL=16
      real tn = fullUpdateTime;

      precomputeStressFromQInterpolated(NorStressGP, XYStressGP, XZStressGP,
          QInterpolatedPlus[face], QInterpolatedMinus[face],  face);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        tn=tn + DeltaT[iTimeGP];

        calcSlipRateAndTraction(TractionGP_XY[iTimeGP], TractionGP_XZ[iTimeGP], NorStressGP[iTimeGP], XYStressGP[iTimeGP], XZStressGP[iTimeGP],
            LocSlipRate, DeltaT[iTimeGP], face);

        calcStateVariableAndFrictionFunc(tmpSlip, LocSlipRate, tn, fullUpdateTime, DeltaT[iTimeGP], face);
      }//End of iTimeGP-Loop

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(LocSlipRate,fullUpdateTime, face);

      //output time when shear stress is equal to the dynamic stress after rupture arrived
      //currently only for linear slip weakening
      outputDynamicStress(fullUpdateTime, face);

      calcPeakSlipRate(LocSlipRate, face);

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      calcAverageSlip(tmpSlip, face);

      postcomputeImposedStateFromNewStress(
          QInterpolatedPlus[face], QInterpolatedMinus[face],
          NorStressGP, TractionGP_XY, TractionGP_XZ,
          timeWeights, face);

      /*
      //debugging
      for(int i = 0; i < tensor::QInterpolated::size(); i++){
        assert( !std::isnan(imposedStatePlus[face][i]) );
      }
      //*/
//      assert(face != 4);

    }//End of Loop over Faces
  }//End of Function evaluate

};//End of Class



class seissol::dr::fr_law::FL_16 : public seissol::dr::fr_law::FL_2 {
public:

  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real tn, real fullUpdateTime,  unsigned int iBndGP,  unsigned int face) override{
    real f2 = 0.0;

    if (lts_t_0[face] == 0) {
      if (tn >= forced_rupture_time[face][iBndGP] ) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max(0.0, std::min( 1.0, (fullUpdateTime - forced_rupture_time[face][iBndGP] ) / lts_t_0[face]));
    }
    stateVariablePsi[iBndGP] = std::max(stateVariablePsi[iBndGP], f2);
  }
};


class seissol::dr::fr_law::FL_3 : public seissol::dr::fr_law::Base {
protected:
  virtual real calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0) {
    return SV0*exp(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-exp(-tmp*time_inc/RS_sl0));
  }

public:
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER],
                        real DeltaT[CONVERGENCE_ORDER]) override {

    seissol::initializers::DR_FL_3 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 *>(dynRup);

    seissol::model::IsotropicWaveSpeeds *waveSpeedsPlus                           = layerData.var(ConcreteLts->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds *waveSpeedsMinus                          = layerData.var(ConcreteLts->waveSpeedsMinus);
    real                    (*initialStressInFaultCS)[numOfPointsPadded][6]       = layerData.var(ConcreteLts->initialStressInFaultCS);
    real                    (*cohesion)[numOfPointsPadded]                        = layerData.var(ConcreteLts->cohesion);
    real*                   RS_f0                                                 = layerData.var(ConcreteLts->RS_f0);
    real*                   RS_a                                                  = layerData.var(ConcreteLts->RS_a);
    real*                   RS_b                                                  = layerData.var(ConcreteLts->RS_b);
    real*                   RS_sl0                                                = layerData.var(ConcreteLts->RS_sl0);
    real*                   RS_sr0                                                = layerData.var(ConcreteLts->RS_sr0);

    real                    (*mu)[numOfPointsPadded]                              = layerData.var(ConcreteLts->mu);
    real                    (*slip)[numOfPointsPadded]                            = layerData.var(ConcreteLts->slip);
    real                    (*slip1)[numOfPointsPadded]                           = layerData.var(ConcreteLts->slip1);
    real                    (*slip2)[numOfPointsPadded]                           = layerData.var(ConcreteLts->slip2);
    real                    (*slipRate1)[numOfPointsPadded]                       = layerData.var(ConcreteLts->slipRate1);
    real                    (*slipRate2)[numOfPointsPadded]                       = layerData.var(ConcreteLts->slipRate2);
    real                    (*rupture_time)[numOfPointsPadded]                    = layerData.var(ConcreteLts->rupture_time);
    bool                    (*RF)[numOfPointsPadded]                              = layerData.var(ConcreteLts->RF);
    real                    (*peakSR)[numOfPointsPadded]                          = layerData.var(ConcreteLts->peakSR);
    real                    (*StateVar)[numOfPointsPadded]                        = layerData.var(ConcreteLts->StateVar);

    real                    (*tracXY)[numOfPointsPadded]                          = layerData.var(ConcreteLts->tracXY);
    real                    (*tracXZ)[numOfPointsPadded]                          = layerData.var(ConcreteLts->tracXZ);
    real                    (*imposedStatePlus)[tensor::QInterpolated::size()]    = layerData.var(ConcreteLts->imposedStatePlus);
    real                    (*imposedStateMinus)[tensor::QInterpolated::size()]   = layerData.var(ConcreteLts->imposedStateMinus);

    //loop parameter are fixed, not variable??
    unsigned int nSRupdates, nSVupdates;
    nSRupdates = 5;
    nSVupdates = 2;



#ifdef _OPENMP
#pragma omp parallel for schedule(static) //private(QInterpolatedPlus,QInterpolatedMinus)
#endif
    for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {

      //TODO: merge TractionGP_XY and tracXY in one variable
      real TractionGP_XY[CONVERGENCE_ORDER][numOfPointsPadded] = {{}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real TractionGP_XZ[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};// OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
      real XYStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
      real XZStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};

      precomputeStressFromQInterpolated(NorStressGP, XYStressGP, XZStressGP,
                                        QInterpolatedPlus[face], QInterpolatedMinus[face], face);

      real LocSlip, LocSlip1, LocSlip2, LocSR1, LocSR2, LocSV, LocCohesion, P_0, LocP, time_inc, P, TotalShearStressYZ, SV0, tmp,tmp2, SlipRateGuess, NR, dNR, LocMu;
      real LocTracXY, LocTracXZ;
      std::array<real, numOfPointsPadded> LocSlipRate;

      for(int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

        LocSlip = slip[face][iBndGP]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
        LocSlip1 = slip1[face][iBndGP]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
        LocSlip2 = slip2[face][iBndGP]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
        LocSR1 = slipRate1[face][iBndGP]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSR2 = slipRate2[face][iBndGP]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSV = StateVar[face][iBndGP];     //DISC%DynRup%StateVar(iBndGP,iFace)
        LocCohesion = cohesion[face][iBndGP]; //DISC%DynRup%cohesion(iBndGP,iFace)          // !< cohesion at given fault node  (should be negative since negative normal stress is compression)
        P_0 = initialStressInFaultCS[face][iBndGP][0]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];

        for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {
          LocP = NorStressGP[iTimeGP][iBndGP];
          time_inc = DeltaT[iTimeGP];

          //SignSR1   = SIGN(1.0,LocSR1)                    ! Gets the sign of the slip rate
          //SignSR2   = SIGN(1.0,LocSR2)                    ! Gets the sign of the slip rate

          // load traction and normal stress
          P = LocP + P_0;

          TotalShearStressYZ = std::sqrt(
              seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][3] + XYStressGP[iTimeGP][iBndGP], 2) +
              seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][5] + XZStressGP[iTimeGP][iBndGP], 2));

          // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996) //TODO: look up
          // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

          SV0 = LocSV;    // Careful, the SV must always be corrected using SV0 and not LocSV!

          // The following process is adapted from that described by Kaneko et al. (2008) TODO: look up

          LocSlipRate[iBndGP]      = std::sqrt(seissol::dr::aux::power(LocSR1,2) + seissol::dr::aux::power(LocSR2,2));
          tmp        = abs( LocSlipRate[iBndGP]);

          for(int j = 0; j < nSVupdates; j++){ //!This loop corrects SV values
            LocSlipRate[iBndGP]=abs( LocSlipRate[iBndGP]);

            // FL= 3 aging law and FL=4 slip law
            LocSV = calcStateVariableHook( SV0,  tmp,  time_inc,  RS_sl0[face]);

            // Newton-Raphson algorithm to determine the value of the slip rate.
            // We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
            //  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
            // In our case we equalize the values of the traction for two equations:
            //             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
            //             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
            //               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta and Rice (2003))

            SlipRateGuess = LocSlipRate[iBndGP];   // SRtest: We use as first guess the SR value of the previous time step

            for(int i = 0; i < nSRupdates; i++){   //!This loop corrects SR values
              tmp          = 0.5/RS_sr0[face]* exp( (RS_f0[face]+RS_b[face]*log(RS_sr0[face]*LocSV/RS_sl0[face]) ) /RS_a[face]);
              tmp2         = tmp * SlipRateGuess;
              NR           = -(1.0/waveSpeedsPlus->sWaveVelocity/waveSpeedsPlus->density+1.0/waveSpeedsMinus->sWaveVelocity/waveSpeedsMinus->density) *
                             (abs(P)*RS_a[face]*log(tmp2+sqrt(seissol::dr::aux::power(tmp2,2)+1.0))-TotalShearStressYZ)-SlipRateGuess;    //!TODO: author before me: not sure if ShTest=TotalShearStressYZ should be + or -...
              dNR          = -(1.0/waveSpeedsPlus->sWaveVelocity/waveSpeedsPlus->density+1.0/waveSpeedsMinus->sWaveVelocity/waveSpeedsMinus->density) *
                             (abs(P)*RS_a[face]/sqrt(1+pow(tmp2,2))*tmp)-1.0;
              SlipRateGuess = abs(SlipRateGuess-NR/dNR);             // no ABS needed around NR/dNR at least for aging law
            }   // End
            tmp=0.5*( LocSlipRate[iBndGP]+abs(SlipRateGuess));  //! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
            LocSlipRate[iBndGP]=abs(SlipRateGuess);
          }   // End SV-Loop

          // FL= 3 aging law and FL=4 slip law
          LocSV= calcStateVariableHook( SV0,  tmp,  time_inc,  RS_sl0[face]);

          //TODO: reused calc from above -> simplify
          tmp  = 0.5 * ( LocSlipRate[iBndGP])/RS_sr0[face] * exp((RS_f0[face] + RS_b[face]*log(RS_sr0[face]*LocSV/RS_sl0[face])) / RS_a[face]);

          LocMu    = RS_a[face] * log(tmp + sqrt(seissol::dr::aux::power(tmp,2) + 1.0));

          // 2D:
          // LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
          // LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results
          // update stress change
          LocTracXY = -((initialStressInFaultCS[face][iBndGP][3] + XYStressGP[iBndGP][iTimeGP])/TotalShearStressYZ)*(LocMu*P+abs(LocCohesion));
          LocTracXZ = -((initialStressInFaultCS[face][iBndGP][5] + XZStressGP[iBndGP][iTimeGP])/TotalShearStressYZ)*(LocMu*P+abs(LocCohesion));
          LocTracXY = LocTracXY - initialStressInFaultCS[face][iBndGP][3];
          LocTracXZ = LocTracXZ - initialStressInFaultCS[face][iBndGP][5];

          // Compute slip
          LocSlip   = LocSlip  + ( LocSlipRate[iBndGP])*time_inc; // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

          //Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
          LocSR1     = -(1.0/(waveSpeedsPlus->sWaveVelocity*waveSpeedsPlus->density)+1.0/(waveSpeedsMinus->sWaveVelocity*waveSpeedsMinus->density))*(LocTracXY-XYStressGP[iTimeGP][iBndGP]);
          LocSR2     = -(1.0/(waveSpeedsPlus->sWaveVelocity*waveSpeedsPlus->density)+1.0/(waveSpeedsMinus->sWaveVelocity*waveSpeedsMinus->density))*(LocTracXZ-XZStressGP[iTimeGP][iBndGP]);

          LocSlip1   = LocSlip1  + (LocSR1)*time_inc;
          LocSlip2   = LocSlip2  + (LocSR2)*time_inc;

          //LocSR1     = SignSR1*ABS(LocSR1)
          //LocSR2     = SignSR2*ABS(LocSR2)

          //Save traction for flux computation
          TractionGP_XY[iTimeGP][iBndGP] = LocTracXY;
          TractionGP_XZ[iTimeGP][iBndGP] = LocTracXZ;
        }//End of iTimeGP- loop

        mu[face][iBndGP]       = LocMu;
        slipRate1[face][iBndGP]  = LocSR1;
        slipRate2[face][iBndGP]  = LocSR2;
        slip[face][iBndGP]       = LocSlip;
        slip1[face][iBndGP]     = LocSlip1;
        slip2[face][iBndGP]      = LocSlip2;
        StateVar[face][iBndGP]   = LocSV;
        tracXY[face][iBndGP] = LocTracXY;
        tracXZ[face][iBndGP] = LocTracXZ;

      }//End of iBndGP-loop

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(LocSlipRate, fullUpdateTime, face);

      calcPeakSlipRate(LocSlipRate, face);


      postcomputeImposedStateFromNewStress(QInterpolatedPlus[face], QInterpolatedMinus[face],
                                           NorStressGP, TractionGP_XY, TractionGP_XZ,
                                           timeWeights, face);
    } //end face-loop
  } //end evaluate function
};

class seissol::dr::fr_law::FL_4 : public seissol::dr::fr_law::FL_3 {
public:

  virtual real calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0) override {
    return RS_sl0/tmp*seissol::dr::aux::power(tmp*SV0/RS_sl0, exp(-tmp*time_inc/RS_sl0));
  }

};


class seissol::dr::fr_law::FL_33 : public seissol::dr::fr_law::Base {
public:
    virtual void evaluate(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup,
                          real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real fullUpdateTime,
                          real timeWeights[CONVERGENCE_ORDER],
                          real DeltaT[CONVERGENCE_ORDER]) override {
        seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
        std::cout << "computing DR for FL_33 (not implemented)\n";
    }
};


#endif //SEISSOL_DR_FRICTION_LAW_H

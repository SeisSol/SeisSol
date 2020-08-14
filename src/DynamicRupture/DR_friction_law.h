//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_FRICTION_LAW_H
#define SEISSOL_DR_FRICTION_LAW_H

#include <c++/8.3.0/iostream>
#include "DR_LTS_Base.h"
#include "DR_math.h"


namespace seissol {
  namespace dr {
    namespace fr_law {
      class Base;
      class FL_2;
      class FL_16;
      class FL_17;
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
  //TODO: implement padded calculation
  //static constexpr int numOfPointsPadded = numberOfPoints;
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];

  /*
   * output:
   * NorStressGP, XYStressGP, XZStressGP
   *
   * input:
   * QInterpolatedPlus, QInterpolatedMinus, eta_p, Zp, Zp_neig, eta_s, Zs, Zs_neig
   */
  void precomputeStressFromQInterpolated(
    real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
    real XYStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
    real XZStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
    real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
    real eta_p, real Zp, real Zp_neig,
    real eta_s, real Zs, real Zs_neig
    ){
    for(int j = 0; j < CONVERGENCE_ORDER; j++){
      auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
      auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      //TODO: does QInterpolatedMinusView work with padded access?
      for(int i = 0; i < numberOfPoints; i++){
        //Carsten Uphoff Thesis: EQ.: 4.53
        NorStressGP[j][i] = eta_p * (QInterpolatedMinusView(i,6) - QInterpolatedPlusView(i,6) + QInterpolatedPlusView(i,0) / Zp + QInterpolatedMinusView(i,0) / Zp_neig);
        XYStressGP[j][i]  = eta_s * (QInterpolatedMinusView(i,7) - QInterpolatedPlusView(i,7) + QInterpolatedPlusView(i,3) / Zs + QInterpolatedMinusView(i,3) / Zs_neig);
        XZStressGP[j][i] = eta_s * (QInterpolatedMinusView(i,8) - QInterpolatedPlusView(i,8) + QInterpolatedPlusView(i,5) / Zs + QInterpolatedMinusView(i,5) / Zs_neig);
      }
    }
    //TODO: is this assert really needed?
    static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],"Different number of quadrature points?");
  }//End of precompute Function

  void postcomputeImposedStateFromNewStress(
      real imposedStatePlus[tensor::QInterpolated::size()],
      real imposedStateMinus[tensor::QInterpolated::size()],
      real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded],
      real TractionGP_XY[CONVERGENCE_ORDER][numOfPointsPadded],
      real TractionGP_XZ[CONVERGENCE_ORDER][numOfPointsPadded],
      real timeWeights[CONVERGENCE_ORDER],
      real Zp, real Zp_neig,
      real Zs, real Zs_neig
      ){
    auto imposedStatePlusView = init::QInterpolated::view::create(imposedStatePlus);
    auto imposedStateMinusView = init::QInterpolated::view::create(imposedStateMinus);
    //initialize to 0
    imposedStateMinusView.setZero();
    imposedStatePlusView.setZero();

    for (int j = 0; j < CONVERGENCE_ORDER; j++) {
      auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
      auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      for (int i = 0; i < numberOfPoints; i++) {
        imposedStateMinusView(i, 0) += timeWeights[j] * NorStressGP[j][i];
        imposedStateMinusView(i, 3) += timeWeights[j] * TractionGP_XY[j][i];
        imposedStateMinusView(i, 5) += timeWeights[j] * TractionGP_XZ[j][i];
        imposedStateMinusView(i, 6) += timeWeights[j] * (QInterpolatedMinusView(i, 6) -  (NorStressGP[j][i] - QInterpolatedMinusView(i, 0))/  Zp_neig);
        imposedStateMinusView(i, 7) += timeWeights[j] * (QInterpolatedMinusView(i, 7) -  (TractionGP_XY[j][i] - QInterpolatedMinusView(i, 3))/  Zs_neig);
        imposedStateMinusView(i, 8) += timeWeights[j] * (QInterpolatedMinusView(i, 8) -  (TractionGP_XZ[j][i] - QInterpolatedMinusView(i, 5))/  Zs_neig);

        imposedStatePlusView(i, 0) += timeWeights[j] * NorStressGP[j][i];
        imposedStatePlusView(i, 3) += timeWeights[j] * TractionGP_XY[j][i];
        imposedStatePlusView(i, 5) += timeWeights[j] * TractionGP_XZ[j][i];
        imposedStatePlusView(i, 6) += timeWeights[j] * (QInterpolatedPlusView(i, 6) +  (NorStressGP[j][i] - QInterpolatedPlusView(i, 0)) /  Zp);
        imposedStatePlusView(i, 7) += timeWeights[j] * (QInterpolatedPlusView(i, 7) +  (TractionGP_XY[j][i] - QInterpolatedPlusView(i, 3)) / Zs);
        imposedStatePlusView(i, 8) += timeWeights[j] * (QInterpolatedPlusView(i, 8) +  (TractionGP_XZ[j][i] - QInterpolatedPlusView(i, 5)) / Zs);
      } //End numberOfPoints-loop
    } //End CONVERGENCE_ORDER-loop

  }

public:
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
public:
  //parameter for insta_healing
  //TODO: make this parameter better accessible?
  real u_0  = 10e-14; //slip rate is considered as being zero for instaneous healing

  //Hook for FL_16
  //TODO: remove this hook - is only needed for evaluate2
  virtual void hook(int iBndGP, real &f1, real tn, real t_0, real fullUpdateTime, real forced_rupture_time[init::QInterpolated::Stop[0]] ) {
      //!Do nothing
  }

  //Hook for FL_16
  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real tn, real t_0, real fullUpdateTime, real forced_rupture_time[init::QInterpolated::Stop[0]] ) {
    //!Do nothing
  }

  //fault strength (Uphoff eq 2.44)
  void calcFaultStrength(std::array<real, numOfPointsPadded> &Strength, real initialStressInFaultCS[numOfPointsPadded][6], real NorStressGP[numOfPointsPadded], real cohesion[numOfPointsPadded], real mu[numOfPointsPadded]){
    std::array<real, numOfPointsPadded> P;
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      P[iBndGP] = initialStressInFaultCS[iBndGP][0] + NorStressGP[iBndGP];
      Strength[iBndGP] = cohesion[iBndGP] - mu[iBndGP] * std::min(P[iBndGP], 0.0);
    }
  }

  void calcTotalShearStressYZ(std::array<real, numOfPointsPadded> &TotalShearStressYZ, real initialStressInFaultCS[numOfPointsPadded][6],
      real XYStressGP[numOfPointsPadded],
      real XZStressGP[numOfPointsPadded]){

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      TotalShearStressYZ[iBndGP] = std::sqrt(
          seissol::dr::aux::power(initialStressInFaultCS[iBndGP][3] + XYStressGP[iBndGP], 2) +
          seissol::dr::aux::power(initialStressInFaultCS[iBndGP][5] + XZStressGP[iBndGP], 2));
    }
  }


  void calcSlipRate(
      real slipRate1[numOfPointsPadded],
      real slipRate2[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real initialStressInFaultCS[numOfPointsPadded][6],
      real XYStressGP[numOfPointsPadded],
      real XZStressGP[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &Strength,
      real eta_s
      ){
    std::array<real, numOfPointsPadded> TotalShearStressYZ;
    std::array<real, numOfPointsPadded> LocSR1;
    std::array<real, numOfPointsPadded> LocSR2;

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      calcTotalShearStressYZ(TotalShearStressYZ, initialStressInFaultCS, XYStressGP, XZStressGP);

      LocSlipRate[iBndGP] = std::max(0.0, (TotalShearStressYZ[iBndGP] - Strength[iBndGP]) / eta_s);

      LocSR1[iBndGP] = LocSlipRate[iBndGP] * (initialStressInFaultCS[iBndGP][3] + XYStressGP[iBndGP]) /
                       (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));
      LocSR2[iBndGP] = LocSlipRate[iBndGP] * (initialStressInFaultCS[iBndGP][5] + XZStressGP[iBndGP]) /
                       (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));

      slipRate1[iBndGP] = LocSR1[iBndGP];
      slipRate2[iBndGP] = LocSR2[iBndGP];
    }
  }

  void calcTraction(
      real TractionGP_XY[numOfPointsPadded],
      real TractionGP_XZ[numOfPointsPadded],
      real tracXY[numOfPointsPadded],
      real tracXZ[numOfPointsPadded],
      real XYStressGP[numOfPointsPadded],
      real XZStressGP[numOfPointsPadded],
      real slipRate1[numOfPointsPadded],
      real slipRate2[numOfPointsPadded],
      real eta_s
      ){

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      TractionGP_XY[iBndGP] = XYStressGP[iBndGP] - eta_s * slipRate1[iBndGP];
      TractionGP_XZ[iBndGP] = XZStressGP[iBndGP] - eta_s * slipRate2[iBndGP];
      tracXY[iBndGP] = TractionGP_XY[iBndGP];
      tracXZ[iBndGP] = TractionGP_XY[iBndGP];
    }
  }

  void updateDirectionalSlip(
      real *slip1,
      real *slip2,
      real *slipRate1,
      real *slipRate2,
      real DeltaT
      ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      //Update slip
      slip1[iBndGP] = slip1[iBndGP] + slipRate1[iBndGP] * DeltaT;
      slip2[iBndGP] = slip2[iBndGP] + slipRate2[iBndGP] * DeltaT;
    }
  }

/*
 * tmp slip is used to calculate the averaged_Slip
 */
  void integrateSliprateToGetSlip(
      real slip[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &tmpSlip,
      //yateto::DenseTensorView<2,double,unsigned> resampleMatrixView,
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real DeltaT
  ){
    std::array<real, numOfPointsPadded> resampledSlipRate{0};
    auto resampleMatrixView = init::resample::view::create(const_cast<double *>(init::resample::Values));


    for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
      //matmul[iBndGP] = 0;

      //TODO: does not work with padded Points bc of resampleMatrix is not padded
      for (int j = 0; j < numberOfPoints; j++) {
        //Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees of freedom
        //resampleMatrix first projects LocSR on the two-dimensional basis on the reference triangle with
        //degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial at the quadrature points
        resampledSlipRate[iBndGP] += resampleMatrixView(iBndGP, j) * LocSlipRate[j];
      }
      slip[iBndGP] = slip[iBndGP] + resampledSlipRate[iBndGP] * DeltaT;
      tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSlipRate[iBndGP] * DeltaT;
    }
  }


  /*
   * the state variable is devided by d_C to simplfy the next equations
   */
  void calcStateVariablePsi(
      std::array<real, numOfPointsPadded> &stateVariablePsi,
      real slip[numOfPointsPadded],
      real d_c[numOfPointsPadded]
      ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
      //f1 = state variable
      //TODO: remove this calc and add it to Friction function calculation
      stateVariablePsi[iBndGP] = std::min(std::abs(slip[iBndGP]) / d_c[iBndGP], 1.0);
    }
  }

  /*
   * Carsten Thesis: Eq. 2.45
   * calculate updated mu -> friction law
   *
   */
  void evaluateFrictionFunction(
      real mu[numOfPointsPadded],
      real mu_S[numOfPointsPadded],
      real mu_D[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &stateVariablePsi
      ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      mu[iBndGP] = mu_S[iBndGP] - (mu_S[iBndGP] - mu_D[iBndGP]) * stateVariablePsi[iBndGP];
    }
  }

  void instaHealingResetMuAndSlip(
      real mu[numOfPointsPadded],
      real slip[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real mu_S[numOfPointsPadded]
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (LocSlipRate[iBndGP] < u_0) {
        mu[iBndGP] = mu_S[iBndGP];
        slip[iBndGP] = 0.0;
      }
    }
  }


  // output rupture front
  // outside of iTimeGP loop in order to safe an 'if' in a loop
  // this way, no subtimestep resolution possible
  void outputRuptureFront(
      bool RF[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &LocSlipRate,
      real rupture_time[numOfPointsPadded],
      real fullUpdateTime
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (RF[iBndGP] && LocSlipRate[iBndGP] > 0.001) {
        rupture_time[iBndGP] = fullUpdateTime;
        RF[iBndGP] = false;
      }
    }
  }

  //output time when shear stress is equal to the dynamic stress after rupture arrived
  //currently only for linear slip weakening
  void outputDynamicStress(
      bool DS[numOfPointsPadded],
      real dynStress_time[numOfPointsPadded],
      real rupture_time[numOfPointsPadded],
      real slip[numOfPointsPadded],
      real d_c[numOfPointsPadded],
      real fullUpdateTime
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

      if (rupture_time[iBndGP] > 0.0 &&
          rupture_time[iBndGP] <= fullUpdateTime &&
          DS[iBndGP] &&
          std::abs(slip[iBndGP]) >= d_c[iBndGP]) {
        dynStress_time[iBndGP] = fullUpdateTime;
        DS[iBndGP] = false;
      }
    }
  }


  void calcPeakSlipRate(
      real peakSR[numOfPointsPadded],
      std::array<real, numOfPointsPadded> &LocSlipRate){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (LocSlipRate[iBndGP] > peakSR[iBndGP]) {
        peakSR[iBndGP] = LocSlipRate[iBndGP];
      }
    }
  }


  //---compute and store slip to determine the magnitude of an earthquake ---
  //    to this end, here the slip is computed and averaged per element
  //    in calc_seissol.f90 this value will be multiplied by the element surface
  //    and an output happened once at the end of the simulation
  void calcAverageSlip(
      real &averaged_Slip,
      bool magnitude_out,
      std::array<real, numOfPointsPadded> &tmpSlip
  ){
    real sum_tmpSlip = 0;
    if (magnitude_out) {
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++)
        sum_tmpSlip += tmpSlip[iBndGP];
      averaged_Slip = averaged_Slip + sum_tmpSlip / numberOfPoints;
    }
  }


  virtual void evaluate(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER],
                        real DeltaT[CONVERGENCE_ORDER]) override {
    //TODO: change later to const_cast
    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);

    real                                (*imposedStatePlus)[tensor::QInterpolated::size()] = layerData.var(ConcreteLts->imposedStatePlus);
    real                                (*imposedStateMinus)[tensor::QInterpolated::size()] = layerData.var(ConcreteLts->imposedStateMinus);
    seissol::model::IsotropicWaveSpeeds *waveSpeedsPlus = layerData.var(ConcreteLts->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds *waveSpeedsMinus = layerData.var(ConcreteLts->waveSpeedsMinus);
    real                    (*initialStressInFaultCS)[init::QInterpolated::Stop[0]][6] = layerData.var(ConcreteLts->initialStressInFaultCS);
    real                    (*cohesion)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->cohesion);
    real                    (*mu)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->mu);
    real                    (*slip)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->slip);
    real                    (*slip1)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->slip1);
    real                    (*slip2)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->slip2);
    real                    (*d_c)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->d_c);
    real                    (*mu_S)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->mu_S);
    real                    (*mu_D)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->mu_D);
    bool                    (*inst_healing) = layerData.var(ConcreteLts->inst_healing);
    real                    (*rupture_time)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->rupture_time);
    bool                    (*RF)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->RF);
    bool                    (*DS)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->DS);
    real                    (*peakSR)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->peakSR);
    real                    (*dynStress_time)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->dynStress_time);
    real                    (*tracXY)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->tracXY);
    real                    (*tracXZ)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->tracXZ);
    real                    (*slipRate1)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->slipRate1);
    real                    (*slipRate2)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->slipRate2);
    bool                    (*magnitude_out) = layerData.var(ConcreteLts->magnitude_out);
    real                    (*averaged_Slip) = layerData.var(ConcreteLts->averaged_Slip);

    //only for FL16:
    real                    (*forced_rupture_time)[init::QInterpolated::Stop[0]] = layerData.var(ConcreteLts->forced_rupture_time);
    real                    *lts_t_0 = layerData.var(ConcreteLts->t_0);

    auto resampleMatrixView = init::resample::view::create(const_cast<double *>(init::resample::Values));

    //TODO: for anisotropic case it must be face dependent?
    //calculate Impedances Z
    real Zp, Zs, Zp_neig, Zs_neig, eta_p, eta_s;
    Zp = (waveSpeedsPlus->density * waveSpeedsPlus->pWaveVelocity);
    Zp_neig = (waveSpeedsMinus->density * waveSpeedsMinus->pWaveVelocity);
    Zs = (waveSpeedsPlus->density * waveSpeedsPlus->sWaveVelocity);
    Zs_neig = (waveSpeedsMinus->density * waveSpeedsMinus->sWaveVelocity);

    eta_p = 1.0 / (1.0 / Zp + 1.0 / Zp_neig);
    eta_s = 1.0 / (1.0 / Zs + 1.0 / Zs_neig);
    //a bit more inaccurate? (after my testings)
    //eta_p = Zp * Zp_neig / (Zp + Zp_neig);
    //eta_s = Zs * Zs_neig / (Zs + Zs_neig);


    #ifdef _OPENMP
    //TODO: do parallel
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
          QInterpolatedPlus[face], QInterpolatedMinus[face],  eta_p, Zp, Zp_neig, eta_s, Zs, Zs_neig);

      //declare local variables
      std::array<real, numOfPointsPadded> tmpSlip{0};
      std::array<real, numOfPointsPadded> Strength;

      std::array<real, numOfPointsPadded> LocSlipRate;
      std::array<real, numOfPointsPadded> stateVariablePsi;
      //tn only needed for FL=16
      real tn = fullUpdateTime;

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        //fault strength (Uphoff eq 2.44)
        calcFaultStrength(Strength, initialStressInFaultCS[face], NorStressGP[iTimeGP], cohesion[face], mu[face]);

        calcSlipRate(slipRate1[face], slipRate2[face], LocSlipRate, initialStressInFaultCS[face],
            XYStressGP[iTimeGP], XZStressGP[iTimeGP], Strength, eta_s);

        calcTraction(TractionGP_XY[iTimeGP], TractionGP_XZ[iTimeGP], tracXY[face],tracXZ[face], XYStressGP[iTimeGP], XZStressGP[iTimeGP],
            slipRate1[face], slipRate2[face], eta_s);

        updateDirectionalSlip(slip1[face], slip2[face], slipRate1[face], slipRate2[face], DeltaT[iTimeGP]);

        integrateSliprateToGetSlip(slip[face], tmpSlip, LocSlipRate, DeltaT[iTimeGP]);

        calcStateVariablePsi(stateVariablePsi, slip[face], d_c[face]);

        //hook for FL_16: may calculate a different value for the state variable Psi
        hookCalcStateVariable(stateVariablePsi, tn, lts_t_0[face], fullUpdateTime, forced_rupture_time[face]); //for FL_16

        //Carsten Thesis: Eq. 2.45
        evaluateFrictionFunction(mu[face], mu_S[face], mu_D[face],stateVariablePsi);

        //instantaneous healing option
        if (inst_healing[face] == true) {
          instaHealingResetMuAndSlip( mu[face], slip[face], LocSlipRate, mu_S[face]);
        }
      }//End of iTimeGP-Loop

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(RF[face],LocSlipRate, rupture_time[face],fullUpdateTime);

      //output time when shear stress is equal to the dynamic stress after rupture arrived
      //currently only for linear slip weakening
      outputDynamicStress(DS[face], dynStress_time[face], rupture_time[face], slip[face], d_c[face], fullUpdateTime);

      calcPeakSlipRate(peakSR[face], LocSlipRate);

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      calcAverageSlip(averaged_Slip[face],magnitude_out[face],tmpSlip);

      postcomputeImposedStateFromNewStress(imposedStatePlus[face], imposedStateMinus[face],
          QInterpolatedPlus[face], QInterpolatedMinus[face],
          NorStressGP, TractionGP_XY, TractionGP_XZ,
          timeWeights, Zp, Zp_neig, Zs, Zs_neig);

      for(int i = 0; i < tensor::QInterpolated::size(); i++){
        assert( !std::isnan(imposedStatePlus[face][i]) );
      }

//      assert(face != 4);
    }//End of Loop over Faces
  }//End of Function evaluate

  /*
   * old version without function cals
   * to compare computational time
   */
  virtual void evaluate2(seissol::initializers::Layer&  layerData,
                         seissol::initializers::DynamicRupture *dynRup,
                         real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                         real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                         real fullUpdateTime,
                         real timeWeights[CONVERGENCE_ORDER],
                         real DeltaT[CONVERGENCE_ORDER]) {

    //TODO: change later to const_cast
    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);
    //std::cout << "computing DR for FL_2\n";

    //DRFaceInformation*                    faceInformation                                                  = layerData.var(ConcreteLts->faceInformation);
    real                                (*imposedStatePlus)[tensor::QInterpolated::size()]                  = layerData.var(ConcreteLts->imposedStatePlus);
    real                                (*imposedStateMinus)[tensor::QInterpolated::size()]                 = layerData.var(ConcreteLts->imposedStateMinus);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                                    = layerData.var(ConcreteLts->waveSpeedsPlus);
    seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                                   = layerData.var(ConcreteLts->waveSpeedsMinus);
    real                    (*initialStressInFaultCS)[init::QInterpolated::Stop[0]][6]                      = layerData.var(ConcreteLts->initialStressInFaultCS);
    real                    (*cohesion)[init::QInterpolated::Stop[0]]                                       = layerData.var(ConcreteLts->cohesion);
    real                    (*mu)[init::QInterpolated::Stop[0]]                                             = layerData.var(ConcreteLts->mu);
    real                    (*slip)[init::QInterpolated::Stop[0]]                                           = layerData.var(ConcreteLts->slip);
    real                    (*slip1)[init::QInterpolated::Stop[0]]                                          = layerData.var(ConcreteLts->slip1);
    real                    (*slip2)[init::QInterpolated::Stop[0]]                                          = layerData.var(ConcreteLts->slip2);
    real                    (*d_c)[init::QInterpolated::Stop[0]]                                            = layerData.var(ConcreteLts->d_c);
    real                    (*mu_S)[init::QInterpolated::Stop[0]]                                           = layerData.var(ConcreteLts->mu_S);
    real                    (*mu_D)[init::QInterpolated::Stop[0]]                                           = layerData.var(ConcreteLts->mu_D);
    bool                    (*inst_healing)                                                                 = layerData.var(ConcreteLts->inst_healing);
    real                    (*rupture_time)[init::QInterpolated::Stop[0]]                                   = layerData.var(ConcreteLts->rupture_time);
    bool                    (*RF)[init::QInterpolated::Stop[0]]                                             = layerData.var(ConcreteLts->RF);
    bool                    (*DS)[init::QInterpolated::Stop[0]]                                             = layerData.var(ConcreteLts->DS);
    real                    (*peakSR)[init::QInterpolated::Stop[0]]                                         = layerData.var(ConcreteLts->peakSR);
    real                    (*dynStress_time)[init::QInterpolated::Stop[0]]                                 = layerData.var(ConcreteLts->dynStress_time);
    real                    (*tracXY)[init::QInterpolated::Stop[0]]                                         = layerData.var(ConcreteLts->tracXY);
    real                    (*tracXZ)[init::QInterpolated::Stop[0]]                                         = layerData.var(ConcreteLts->tracXZ);
    real                    (*slipRate1)[init::QInterpolated::Stop[0]]                                      = layerData.var(ConcreteLts->slipRate1);
    real                    (*slipRate2)[init::QInterpolated::Stop[0]]                                      = layerData.var(ConcreteLts->slipRate2);
    bool                    (*magnitude_out)                                                                = layerData.var(ConcreteLts->magnitude_out);
    real                    (*averaged_Slip)                                                                = layerData.var(ConcreteLts->averaged_Slip);

    //only for FL16:
    real                    (*forced_rupture_time)[init::QInterpolated::Stop[0]]                            = layerData.var(ConcreteLts->forced_rupture_time);
    real*                                   lts_t_0                                                         = layerData.var(ConcreteLts->t_0);

    auto resampleMatrixView = init::resample::view::create(const_cast<double *>(init::resample::Values));



    //TODO: for anisotropic case it must be face dependent?
    //calculate Impedances Z
    real Zp, Zs, Zp_neig, Zs_neig, eta_p, eta_s;
    Zp = (waveSpeedsPlus->density * waveSpeedsPlus->pWaveVelocity);
    Zp_neig = (waveSpeedsMinus->density * waveSpeedsMinus->pWaveVelocity);
    Zs = (waveSpeedsPlus->density * waveSpeedsPlus->sWaveVelocity);
    Zs_neig = (waveSpeedsMinus->density * waveSpeedsMinus->sWaveVelocity);

    eta_p = 1.0 / (1.0/Zp + 1.0/Zp_neig);
    eta_s = 1.0 / (1.0/Zs + 1.0/Zs_neig);
    //a bit more inaccurate? (after my testings)
    //eta_p = Zp * Zp_neig / (Zp + Zp_neig);
    //eta_s = Zs * Zs_neig / (Zs + Zs_neig);

#ifdef _OPENMP
//TODO do parallel
#pragma omp parallel for schedule(static) //private(QInterpolatedPlus,QInterpolatedMinus)
#endif
    //TODO: split loop
    for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {

      real TractionGP_XY2[CONVERGENCE_ORDER][numOfPointsPadded] = {{}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real TractionGP_XZ2[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};// OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real NorStressGP2[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
      real XYStressGP2[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
      real XZStressGP2[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};

      precomputeStressFromQInterpolated(NorStressGP2, XYStressGP2, XZStressGP2,
                                        QInterpolatedPlus[face], QInterpolatedMinus[face],  eta_p, Zp, Zp_neig, eta_s, Zs, Zs_neig);
/*
      real TractionGP_XY[numOfPointsPadded][CONVERGENCE_ORDER] = {{}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real TractionGP_XZ[numOfPointsPadded][CONVERGENCE_ORDER] = {{}};// OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
      real NorStressGP[numOfPointsPadded][CONVERGENCE_ORDER] = {{}};
      real XYStressGP[numOfPointsPadded][CONVERGENCE_ORDER]= {{}};
      real XZStressGP[numOfPointsPadded][CONVERGENCE_ORDER]= {{}};

      for(int j = 0; j < CONVERGENCE_ORDER; j++){
        auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[face][j]);
        auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[face][j]);
        //TODO: does QInterpolatedMinusView work with padded access?
        for(int i = 0; i < numOfPointsPadded; i++){
          //Carsten Uphoff Thesis: EQ.: 4.53
          NorStressGP[i][j] = eta_p * (QInterpolatedMinusView(i,6) - QInterpolatedPlusView(i,6) + QInterpolatedPlusView(i,0) / Zp + QInterpolatedMinusView(i,0) / Zp_neig);
          XYStressGP[i][j]  = eta_s * (QInterpolatedMinusView(i,7) - QInterpolatedPlusView(i,7) + QInterpolatedPlusView(i,3) / Zs + QInterpolatedMinusView(i,3) / Zs_neig);
          XZStressGP[i][j] = eta_s * (QInterpolatedMinusView(i,8) - QInterpolatedPlusView(i,8) + QInterpolatedPlusView(i,5) / Zs + QInterpolatedMinusView(i,5) / Zs_neig);
        }
      }
*/

      //TODO: is this assert really needed?
      static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0], "Different number of quadrature points?");

      //required input:

      //declare local variables
      real sum_tmpSlip;
      real stateVariablePsi;
      //real time_inc = 0;

      //Initialized to Zero
      //real tmpSlip[numOfPointsPadded] = {0}; //zero initialized
      //real matmul[numOfPointsPadded] = {0};
      std::array<real, numOfPointsPadded> tmpSlip{0};
      std::array<real, numOfPointsPadded> matmul{0};

      //no initialization required:
      std::array<real, numOfPointsPadded> P;
      std::array<real, numOfPointsPadded> Strength;
      std::array<real, numOfPointsPadded> TotalShearStressYZ;
      std::array<real, numOfPointsPadded> LocSR;
      std::array<real, numOfPointsPadded> LocSR1;
      std::array<real, numOfPointsPadded> LocSR2;
      std::array<real, numOfPointsPadded> LocTracXY;
      std::array<real, numOfPointsPadded> LocTracXZ;

      //tn only needed for FL=16
      real tn = fullUpdateTime;
      //real t_0 = lts_t_0[face]; //= DISC%DynRup%t_0     used for FL16

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
          //time_inc = DeltaT[iTimeGP];
          tn = tn + DeltaT[iTimeGP]; //could be moved to FL_16 hook()

          //P[iBndGP] = initialStressInFaultCS[face][iBndGP][0] + NorStressGP[iBndGP][iTimeGP];
          P[iBndGP] = initialStressInFaultCS[face][iBndGP][0] + NorStressGP2[iTimeGP][iBndGP];

          //fault strength (Uphoff eq 2.44)
          Strength[iBndGP] = cohesion[face][iBndGP] - mu[face][iBndGP] * std::min(P[iBndGP], 0.0);

          //Carsten Uphoff Thesis: Eq. 4.59
          //TotalShearStressYZ[iBndGP] = std::sqrt(
          //    seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][3] + XYStressGP[iBndGP][iTimeGP], 2) +
          //    seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][5] + XZStressGP[iBndGP][iTimeGP], 2));
          TotalShearStressYZ[iBndGP] = std::sqrt(
              seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][3] + XYStressGP2[iTimeGP][iBndGP], 2) +
              seissol::dr::aux::power(initialStressInFaultCS[face][iBndGP][5] + XZStressGP2[iTimeGP][iBndGP], 2));


          LocSR[iBndGP] = std::max(0.0, (TotalShearStressYZ[iBndGP] - Strength[iBndGP]) / eta_s);
          /*
          LocSR1[iBndGP] =
                  LocSR[iBndGP] * (frD.getInitialStressInFaultCS(iBndGP, 3, iFace) + XYStressGP[iBndGP][iTimeGP]) /
                  (Strength[iBndGP] + eta * LocSR[iBndGP]);
          LocSR2[iBndGP] =
                  LocSR[iBndGP] * (frD.getInitialStressInFaultCS(iBndGP, 5, iFace) + XZStressGP[iBndGP][iTimeGP]) /
                  (Strength[iBndGP] + eta * LocSR[iBndGP]);
          */
          //TODO: check alternative faster calc??
          /*
          LocSR1[iBndGP] = LocSR[iBndGP] * (initialStressInFaultCS[face][iBndGP][3] +
                                            XYStressGP[iBndGP][iTimeGP]) /
                           (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));
          LocSR2[iBndGP] = LocSR[iBndGP] * (initialStressInFaultCS[face][iBndGP][5] +
                                            XZStressGP[iBndGP][iTimeGP]) /
                           (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));


          LocTracXY[iBndGP] = XYStressGP[iBndGP][iTimeGP] - eta_s * LocSR1[iBndGP];
          LocTracXZ[iBndGP] = XZStressGP[iBndGP][iTimeGP] - eta_s * LocSR2[iBndGP];
          */

          LocSR1[iBndGP] = LocSR[iBndGP] * (initialStressInFaultCS[face][iBndGP][3] +
                                            XYStressGP2[iTimeGP][iBndGP]) /
                           (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));
          LocSR2[iBndGP] = LocSR[iBndGP] * (initialStressInFaultCS[face][iBndGP][5] +
                                            XZStressGP2[iTimeGP][iBndGP]) /
                           (std::max(TotalShearStressYZ[iBndGP], Strength[iBndGP]));
          LocTracXY[iBndGP] = XYStressGP2[iTimeGP][iBndGP] - eta_s * LocSR1[iBndGP];
          LocTracXZ[iBndGP] = XZStressGP2[iTimeGP][iBndGP] - eta_s * LocSR2[iBndGP];

          //Update slip
          slip1[face][iBndGP] = slip1[face][iBndGP] + LocSR1[iBndGP] * DeltaT[iTimeGP];
          slip2[face][iBndGP] = slip2[face][iBndGP] + LocSR2[iBndGP] * DeltaT[iTimeGP];
        }
        //TODO: does not work with padded Points bc of resampleMatrix is not padded
        for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
          //Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees of freedom
          //resampleMatrix first projects LocSR on the two-dimensional basis on the reference triangle with
          //degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial at the quadrature points
          matmul[iBndGP] = 0;

          //TODO: does not work with padded Points bc of resampleMatrix is not padded
          for (int j = 0; j < numberOfPoints; j++) {
            matmul[iBndGP] += resampleMatrixView(iBndGP, j) * LocSR[j];
          }
          slip[face][iBndGP] = slip[face][iBndGP] + matmul[iBndGP] * DeltaT[iTimeGP];
          tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSR[iBndGP] * DeltaT[iTimeGP];


          //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
          //f1 = state variable
          stateVariablePsi = std::min(std::abs(slip[face][iBndGP]) / d_c[face][iBndGP], 1.0);

          //hook for FL_16: may calculate a different value for the state variable f1
          hook(iBndGP, stateVariablePsi, tn, lts_t_0[face], fullUpdateTime, forced_rupture_time[init::QInterpolated::Stop[0]]); //for FL_16

          //Carsten Thesis: Eq. 2.45
          mu[face][iBndGP] = mu_S[face][iBndGP] - (mu_S[face][iBndGP] - mu_D[face][iBndGP]) * stateVariablePsi;

          //instantaneous healing option
          if (inst_healing[face] == true) {
            if (LocSR[iBndGP] < u_0) {
              mu[face][iBndGP] = mu_S[face][iBndGP];
              slip[face][iBndGP] = 0.0;
            }
          }

          //TODO: why do we need LocTracXY
          TractionGP_XY2[iTimeGP][iBndGP] = LocTracXY[iBndGP];
          TractionGP_XZ2[iTimeGP][iBndGP] = LocTracXZ[iBndGP];

        }
      }     //end iTimeGP loop

      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {

        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if (RF[face][iBndGP] && LocSR[iBndGP] > 0.001) {
          rupture_time[face][iBndGP] = fullUpdateTime;
          RF[face][iBndGP] = false;
        }
        //output time when shear stress is equal to the dynamic stress after rupture arrived
        //currently only for linear slip weakening
        if (rupture_time[face][iBndGP] > 0.0 &&
            rupture_time[face][iBndGP] <= fullUpdateTime &&
            DS[face][iBndGP] &&
            std::abs(slip[face][iBndGP]) >= d_c[face][iBndGP]) {
          dynStress_time[face][iBndGP] = fullUpdateTime;
          DS[face][iBndGP] = false;
        }

        if (LocSR[iBndGP] > peakSR[face][iBndGP]) {
          peakSR[face][iBndGP] = LocSR[iBndGP];
        }

        tracXY[face][iBndGP] = LocTracXY[iBndGP];
        tracXZ[face][iBndGP] = LocTracXZ[iBndGP];

        slipRate1[face][iBndGP] = LocSR1[iBndGP];
        slipRate2[face][iBndGP] = LocSR2[iBndGP];
      } //end i < nBndGP loop

      //---compute and store slip to determine the magnitude of an earthquake ---
      //    to this end, here the slip is computed and averaged per element
      //    in calc_seissol.f90 this value will be multiplied by the element surface
      //    and an output happened once at the end of the simulation
      sum_tmpSlip = 0;
      if (magnitude_out[face]) {
        for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
          sum_tmpSlip += tmpSlip[iBndGP];
        }
        averaged_Slip[face] = averaged_Slip[face] + sum_tmpSlip / numberOfPoints;
      }

      postcomputeImposedStateFromNewStress(imposedStatePlus[face], imposedStateMinus[face],
                                           QInterpolatedPlus[face], QInterpolatedMinus[face],
                                           NorStressGP2, TractionGP_XY2, TractionGP_XZ2,
                                           timeWeights, Zp, Zp_neig, Zs, Zs_neig);
      /*
      auto imposedStatePlusView = init::QInterpolated::view::create(imposedStatePlus[face]);
      auto imposedStateMinusView = init::QInterpolated::view::create(imposedStateMinus[face]);
      //initialize to 0
      imposedStateMinusView.setZero();
      imposedStatePlusView.setZero();

      for (int j = 0; j < CONVERGENCE_ORDER; j++) {
        auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[face][j]);
        auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[face][j]);
        for (int i = 0; i < numberOfPoints; i++) {
          imposedStateMinusView(i, 0) += timeWeights[j] * NorStressGP[i][j];
          imposedStateMinusView(i, 3) += timeWeights[j] * TractionGP_XY[i][j];
          imposedStateMinusView(i, 5) += timeWeights[j] * TractionGP_XZ[i][j];
          imposedStateMinusView(i, 6) += timeWeights[j] * (QInterpolatedMinusView(i, 6) -  (NorStressGP[i][j] - QInterpolatedMinusView(i, 0))/  Zp_neig);
          imposedStateMinusView(i, 7) += timeWeights[j] * (QInterpolatedMinusView(i, 7) -  (TractionGP_XY[i][j] - QInterpolatedMinusView(i, 3))/  Zs_neig);
          imposedStateMinusView(i, 8) += timeWeights[j] * (QInterpolatedMinusView(i, 8) -  (TractionGP_XZ[i][j] - QInterpolatedMinusView(i, 5))/  Zs_neig);

          imposedStatePlusView(i, 0) += timeWeights[j] * NorStressGP[i][j];
          imposedStatePlusView(i, 3) += timeWeights[j] * TractionGP_XY[i][j];
          imposedStatePlusView(i, 5) += timeWeights[j] * TractionGP_XZ[i][j];
          imposedStatePlusView(i, 6) += timeWeights[j] * (QInterpolatedPlusView(i, 6) +  (NorStressGP[i][j] - QInterpolatedPlusView(i, 0)) /  Zp);
          imposedStatePlusView(i, 7) += timeWeights[j] * (QInterpolatedPlusView(i, 7) +  (TractionGP_XY[i][j] - QInterpolatedPlusView(i, 3)) / Zs);
          imposedStatePlusView(i, 8) += timeWeights[j] * (QInterpolatedPlusView(i, 8) +  (TractionGP_XZ[i][j] - QInterpolatedPlusView(i, 5)) / Zs);
        } //End numberOfPoints-loop
      } //End CONVERGENCE_ORDER-loop

      //assert(face != 4);
      */
      /*
      real XXXimposed[468] = {0};

      for(int i = 0; i < tensor::QInterpolated::size(); i++){
          XXXimposed[i] = imposedStateMinus[face][i];
      }
      for (int i = 0; i < numberOfPoints; i++) {
          for(int j = 0; j < 9; j++){
              assert( !std::isnan( imposedStateMinusView(i,j) ) );
          }

      }
      for(int i = 0; i < tensor::QInterpolated::size(); i++){
          assert( !std::isnan(imposedStateMinus[face][i]) );
      }
      */
    }//End of loop over all faces
  }//End of Function evaluate


};//End of Class

class seissol::dr::fr_law::FL_17 : public seissol::dr::fr_law::Base {
public:
    virtual void hook() {}

    virtual void evaluate(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup,
                          real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real fullUpdateTime,
                          real timeWeights[CONVERGENCE_ORDER],
                          real DeltaT[CONVERGENCE_ORDER]) override {
        //std::cout << "computing ";
        hook();
        //std::cout << " DR for FL_16\n";
    }
};

class seissol::dr::fr_law::FL_16 : public seissol::dr::fr_law::FL_2 {
public:
  /*
   * f1 output
   */
  virtual void hook(int iBndGP, real &f1, real tn, real t_0, real fullUpdateTime, real forced_rupture_time[init::QInterpolated::Stop[0]] )  override {
    real f2 = 0.0;

    if (t_0 == 0) {
      if (tn >= forced_rupture_time[iBndGP] ) {
          f2 = 1.0;
      } else {
          f2 = 0.0;
      }
    } else {
      f2 = std::max(0.0, std::min((fullUpdateTime - forced_rupture_time[iBndGP] ) / t_0, 1.0));
    }
    f1 = std::max(f1, f2);
  }

  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real tn, real t_0, real fullUpdateTime, real forced_rupture_time[init::QInterpolated::Stop[0]] ) override{
    real f2 = 0.0;

    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (t_0 == 0) {
        if (tn >= forced_rupture_time[iBndGP] ) {
          f2 = 1.0;
        } else {
          f2 = 0.0;
        }
      } else {
        f2 = std::max(0.0, std::min((fullUpdateTime - forced_rupture_time[iBndGP] ) / t_0, 1.0));
      }
      stateVariablePsi[iBndGP] = std::max(stateVariablePsi[iBndGP], f2);
    }

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
        std::cout << "computing DR for FL_33\n";
    }
};


#endif //SEISSOL_DR_FRICTION_LAW_H

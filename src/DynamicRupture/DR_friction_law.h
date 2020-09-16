//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_FRICTION_LAW_H
#define SEISSOL_DR_FRICTION_LAW_H

#include <c++/8.3.0/iostream>
#include "DR_math.h"
#include <yaml-cpp/yaml.h>


namespace seissol {
  namespace dr {
    namespace fr_law {
      class BaseFrictionSolver;
      class LinearSlipWeakeningSolverFL2;  //linear slip weakening
      class Solver_FL_3; //rate and state aging law
      class Solver_FL_4; //rate and state slip law
      class Solver_FL_16; //Linear slip weakening forced time rapture
      class Solver_FL_33; //ImposedSlipRateOnDRBoundary
      class Solver_FL_103;  //rate and state nuc103
    }
  }
}


class seissol::dr::fr_law::BaseFrictionSolver {

public:
  //TODO: remove (later) only for debugging:
  int numberOfFunctionCalls = 0;

  virtual ~BaseFrictionSolver() {}

  //set the parameters from .par file with yaml to this class attributes.
  void setInputParam(const YAML::Node& Params) {
    using namespace initializers;
    //TODO: maybe allocate dr::DrParameterT in MemoryManager and copy here only the reference
    m_Params.setAllInputParam(Params);
  }


protected:
  static constexpr int numberOfPoints =  tensor::QInterpolated::Shape[0];// DISC%Galerkin%nBndGP
  //TODO: is init::QInterpolated::Start[0] always 0?
  //assert(init::QInterpolated::Start[0] == 0);
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
  //YAML::Node m_InputParam;
  dr::DrParameterT m_Params;
  ImpedancesAndEta*                     impAndEta;
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

  struct FaultStresses{
    //TODO: merge TractionGP_XY and tracXY in one variable
    real TractionGP_XY[CONVERGENCE_ORDER][numOfPointsPadded] = {{}}; // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
    real TractionGP_XZ[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};// OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
    real NorStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
    real XYStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
    real XZStressGP[CONVERGENCE_ORDER][numOfPointsPadded] = {{}};
  };

  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  virtual void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup){
    impAndEta                                     = layerData.var(dynRup->impAndEta);
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
    FaultStresses &faultStresses,
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
        faultStresses.NorStressGP[j][i] = impAndEta[face].eta_p * (QInterpolatedMinusView(i,6) - QInterpolatedPlusView(i,6) + QInterpolatedPlusView(i,0) / impAndEta[face].Zp + QInterpolatedMinusView(i,0) / impAndEta[face].Zp_neig);
        faultStresses.XYStressGP[j][i]  = impAndEta[face].eta_s * (QInterpolatedMinusView(i,7) - QInterpolatedPlusView(i,7) + QInterpolatedPlusView(i,3) / impAndEta[face].Zs + QInterpolatedMinusView(i,3) / impAndEta[face].Zs_neig);
        faultStresses.XZStressGP[j][i] = impAndEta[face].eta_s * (QInterpolatedMinusView(i,8) - QInterpolatedPlusView(i,8) + QInterpolatedPlusView(i,5) / impAndEta[face].Zs + QInterpolatedMinusView(i,5) / impAndEta[face].Zs_neig);
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
      const FaultStresses &faultStresses,
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
        imposedStateMinusView(i, 0) += timeWeights[j] * faultStresses.NorStressGP[j][i];
        imposedStateMinusView(i, 3) += timeWeights[j] * faultStresses.TractionGP_XY[j][i];
        imposedStateMinusView(i, 5) += timeWeights[j] * faultStresses.TractionGP_XZ[j][i];
        imposedStateMinusView(i, 6) += timeWeights[j] * (QInterpolatedMinusView(i, 6) -  (faultStresses.NorStressGP[j][i] - QInterpolatedMinusView(i, 0))/  impAndEta[face].Zp_neig);
        imposedStateMinusView(i, 7) += timeWeights[j] * (QInterpolatedMinusView(i, 7) -  (faultStresses.TractionGP_XY[j][i] - QInterpolatedMinusView(i, 3))/  impAndEta[face].Zs_neig);
        imposedStateMinusView(i, 8) += timeWeights[j] * (QInterpolatedMinusView(i, 8) -  (faultStresses.TractionGP_XZ[j][i] - QInterpolatedMinusView(i, 5))/  impAndEta[face].Zs_neig);

        imposedStatePlusView(i, 0) += timeWeights[j] * faultStresses.NorStressGP[j][i];
        imposedStatePlusView(i, 3) += timeWeights[j] * faultStresses.TractionGP_XY[j][i];
        imposedStatePlusView(i, 5) += timeWeights[j] * faultStresses.TractionGP_XZ[j][i];
        imposedStatePlusView(i, 6) += timeWeights[j] * (QInterpolatedPlusView(i, 6) +  (faultStresses.NorStressGP[j][i] - QInterpolatedPlusView(i, 0)) /  impAndEta[face].Zp);
        imposedStatePlusView(i, 7) += timeWeights[j] * (QInterpolatedPlusView(i, 7) +  (faultStresses.TractionGP_XY[j][i] - QInterpolatedPlusView(i, 3)) / impAndEta[face].Zs);
        imposedStatePlusView(i, 8) += timeWeights[j] * (QInterpolatedPlusView(i, 8) +  (faultStresses.TractionGP_XZ[j][i] - QInterpolatedPlusView(i, 5)) / impAndEta[face].Zs);
      } //End numberOfPoints-loop
    } //End CONVERGENCE_ORDER-loop

  }

  // output rupture front
  // outside of iTimeGP loop in order to safe an 'if' in a loop
  // this way, no subtimestep resolution possible
  void outputRuptureFront(
      real LocSlipRate[seissol::tensor::resamplePar::size()],
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
      real LocSlipRate[seissol::tensor::resamplePar::size()],
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
class seissol::dr::fr_law::LinearSlipWeakeningSolverFL2 : public seissol::dr::fr_law::BaseFrictionSolver {
protected:
  //parameter for insta_healing
  //TODO: make this parameter better accessible?
  real u_0  = 10e-14; //slip rate is considered as being zero for instaneous healing
  real                    (*d_c)[numOfPointsPadded];
  real                    (*mu_S)[numOfPointsPadded];
  real                    (*mu_D)[numOfPointsPadded];
  bool                    (*DS)[numOfPointsPadded];
  real                    (*averaged_Slip);
  real                    (*dynStress_time)[numOfPointsPadded];
  //only for FL16:
  real                    (*forced_rupture_time)[numOfPointsPadded];

  //Hook for Factory_FL_16
  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real tn, real fullUpdateTime, unsigned int iBndGP, unsigned int face  ) {
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
  }

  /*
   *
   */
  void calcSlipRateAndTraction(
      FaultStresses faultStresses,
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
      real tn,
      real fullUpdateTime,
      real DeltaT,
      unsigned int face
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
      slip[face][iBndGP] = slip[face][iBndGP] +  resampledSlipRate[iBndGP] * DeltaT;
      tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSlipRate[iBndGP] * DeltaT;

      //-------------------------------------
      //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
      //actually slip is already the stateVariable for this FL, but to simplify the next equations we divide it here by d_C
      stateVariablePsi[iBndGP] = std::min(std::fabs(slip[face][iBndGP]) / d_c[face][iBndGP], 1.0);

      //-------------------------------------
      //hook for Factory_FL_16: may calculate a different value for the state variable Psi
      hookCalcStateVariable(stateVariablePsi, tn, fullUpdateTime, iBndGP, face); //for Factory_FL_16

      //-------------------------------------
      //Carsten Thesis: Eq. 2.45
      //evaluate friction law: updated mu -> friction law
      mu[face][iBndGP] = mu_S[face][iBndGP] - (mu_S[face][iBndGP] - mu_D[face][iBndGP]) * stateVariablePsi[iBndGP];

      //-------------------------------------
      //instantaneous healing option Reset Mu and Slip
      if (m_Params.IsInstaHealingOn == true) {
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
          std::fabs(slip[face][iBndGP]) >= d_c[face][iBndGP]) {
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
    if (m_Params.IsMagnitudeOutputOn) {
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
      //initialize struct for in/outputs stresses
      FaultStresses faultStresses{};

      //declare local variables
      std::array<real, numOfPointsPadded> tmpSlip{0};
      real LocSlipRate[seissol::tensor::resamplePar::size()];
      dynamicRupture::kernel::resampleParameter resampleKrnl;
      resampleKrnl.resampleM = init::resample::Values;
      //tn only needed for FL=16
      real tn = fullUpdateTime;

      precomputeStressFromQInterpolated(faultStresses,
          QInterpolatedPlus[face], QInterpolatedMinus[face],  face);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        tn=tn + DeltaT[iTimeGP];
        calcSlipRateAndTraction(faultStresses, LocSlipRate, DeltaT[iTimeGP], iTimeGP ,face);

        calcStateVariableAndFrictionFunc(tmpSlip, LocSlipRate, resampleKrnl, tn, fullUpdateTime, DeltaT[iTimeGP], face);
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
          faultStresses, timeWeights, face);
    }//End of Loop over Faces
  }//End of Function evaluate

};//End of Class



class seissol::dr::fr_law::Solver_FL_16 : public seissol::dr::fr_law::LinearSlipWeakeningSolverFL2 {
public:

  virtual void hookCalcStateVariable(std::array<real, numOfPointsPadded> &stateVariablePsi, real tn, real fullUpdateTime,  unsigned int iBndGP,  unsigned int face) override{
    real f2 = 0.0;

    if (m_Params.t_0 == 0) {
      if (tn >= forced_rupture_time[face][iBndGP] ) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max(0.0, std::min( 1.0, (fullUpdateTime - forced_rupture_time[face][iBndGP] ) / m_Params.t_0));
    }
    stateVariablePsi[iBndGP] = std::max(stateVariablePsi[iBndGP], f2);
  }
};



class seissol::dr::fr_law::Solver_FL_33 : public seissol::dr::fr_law::BaseFrictionSolver {
public:
    virtual void evaluate(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup,
                          real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                          real fullUpdateTime,
                          real timeWeights[CONVERGENCE_ORDER],
                          real DeltaT[CONVERGENCE_ORDER]) override {
        seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
        std::cout << "computing DR for Init_FL_33 (not implemented)\n";
    }
};



#endif //SEISSOL_DR_FRICTION_LAW_H

//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_SOLVER_BASE_H
#define SEISSOL_DR_SOLVER_BASE_H

#include <c++/8.3.0/iostream>
#include "DR_math.h"
#include <yaml-cpp/yaml.h>


namespace seissol {
  namespace dr {
    namespace fr_law {
      class BaseFrictionSolver;
      class SolverNoFaultFL0;
      class Solver_FL_33; //ImposedSlipRateOnDRBoundary
      class SolverBluePrint;
    }
  }
}


class seissol::dr::fr_law::BaseFrictionSolver {

public:
  //TODO: remove (later) only for debugging:
  int numberOfFunctionCalls = 0;

  virtual ~BaseFrictionSolver() {}

  //set the parameters from .par file with yaml to this class attributes.
  void setInputParam(dr::DrParameterT *DynRupParameter) {
    m_Params = DynRupParameter;
  }


protected:
  static constexpr int numberOfPoints =  tensor::QInterpolated::Shape[0];// DISC%Galerkin%nBndGP
  //TODO: is init::QInterpolated::Start[0] always 0?
  //assert(init::QInterpolated::Start[0] == 0);
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
  //YAML::Node m_InputParam;
  dr::DrParameterT *m_Params;
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
  real                    (*tracXY)[numOfPointsPadded];
  real                    (*tracXZ)[numOfPointsPadded];
  real                    (*imposedStatePlus)[tensor::QInterpolated::size()];
  real                    (*imposedStateMinus)[tensor::QInterpolated::size()];

  //only for some FLs initialized:
  real  *averaged_Slip;

  struct FaultStresses{
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
    unsigned int ltsFace
    ){

    for(int j = 0; j < CONVERGENCE_ORDER; j++){
      auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
      auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      for(int i = 0; i < numberOfPoints; i++){
        //Carsten Uphoff Thesis: EQ.: 4.53
        faultStresses.NorStressGP[j][i] = impAndEta[ltsFace].eta_p * (QInterpolatedMinusView(i, 6) - QInterpolatedPlusView(i, 6) + QInterpolatedPlusView(i, 0) / impAndEta[ltsFace].Zp + QInterpolatedMinusView(i, 0) / impAndEta[ltsFace].Zp_neig);
        faultStresses.XYStressGP[j][i]  = impAndEta[ltsFace].eta_s * (QInterpolatedMinusView(i, 7) - QInterpolatedPlusView(i, 7) + QInterpolatedPlusView(i, 3) / impAndEta[ltsFace].Zs + QInterpolatedMinusView(i, 3) / impAndEta[ltsFace].Zs_neig);
        faultStresses.XZStressGP[j][i] = impAndEta[ltsFace].eta_s * (QInterpolatedMinusView(i, 8) - QInterpolatedPlusView(i, 8) + QInterpolatedPlusView(i, 5) / impAndEta[ltsFace].Zs + QInterpolatedMinusView(i, 5) / impAndEta[ltsFace].Zs_neig);
      }
    }
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
      unsigned int ltsFace
      ){
    auto imposedStatePlusView = init::QInterpolated::view::create(imposedStatePlus[ltsFace]);
    auto imposedStateMinusView = init::QInterpolated::view::create(imposedStateMinus[ltsFace]);
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
        imposedStateMinusView(i, 6) += timeWeights[j] * (QInterpolatedMinusView(i, 6) -  (faultStresses.NorStressGP[j][i] - QInterpolatedMinusView(i, 0))/  impAndEta[ltsFace].Zp_neig);
        imposedStateMinusView(i, 7) += timeWeights[j] * (QInterpolatedMinusView(i, 7) -  (faultStresses.TractionGP_XY[j][i] - QInterpolatedMinusView(i, 3))/  impAndEta[ltsFace].Zs_neig);
        imposedStateMinusView(i, 8) += timeWeights[j] * (QInterpolatedMinusView(i, 8) -  (faultStresses.TractionGP_XZ[j][i] - QInterpolatedMinusView(i, 5))/  impAndEta[ltsFace].Zs_neig);

        imposedStatePlusView(i, 0) += timeWeights[j] * faultStresses.NorStressGP[j][i];
        imposedStatePlusView(i, 3) += timeWeights[j] * faultStresses.TractionGP_XY[j][i];
        imposedStatePlusView(i, 5) += timeWeights[j] * faultStresses.TractionGP_XZ[j][i];
        imposedStatePlusView(i, 6) += timeWeights[j] * (QInterpolatedPlusView(i, 6) +  (faultStresses.NorStressGP[j][i] - QInterpolatedPlusView(i, 0)) /  impAndEta[ltsFace].Zp);
        imposedStatePlusView(i, 7) += timeWeights[j] * (QInterpolatedPlusView(i, 7) +  (faultStresses.TractionGP_XY[j][i] - QInterpolatedPlusView(i, 3)) / impAndEta[ltsFace].Zs);
        imposedStatePlusView(i, 8) += timeWeights[j] * (QInterpolatedPlusView(i, 8) +  (faultStresses.TractionGP_XZ[j][i] - QInterpolatedPlusView(i, 5)) / impAndEta[ltsFace].Zs);
      } //End numberOfPoints-loop
    } //End CONVERGENCE_ORDER-loop

  }

  /*
  * Function from NucleationFunctions_mod.f90
  */
  double Calc_SmoothStepIncrement(double fullUpdateTime, real dt){
    double Gnuc;
    double prevtime;
    if(fullUpdateTime > 0.0 && fullUpdateTime <= m_Params->t_0){
      Gnuc = Calc_SmoothStep(fullUpdateTime);
      prevtime = fullUpdateTime - dt;
      if(prevtime > 0.0){
        Gnuc = Gnuc - Calc_SmoothStep(prevtime);
      }
    }else{
      Gnuc = 0.0;
    }
    return Gnuc;
  }

  /*
  * Function in NucleationFunctions_mod.f90
  */
  double Calc_SmoothStep(double fullUpdateTime){
    double Gnuc;
    if (fullUpdateTime <= 0){
      Gnuc=0.0;
    }else{
      if (fullUpdateTime < m_Params->t_0){
        Gnuc = std::exp(seissol::dr::aux::power(fullUpdateTime - m_Params->t_0, 2) / (fullUpdateTime * (fullUpdateTime - 2.0 * m_Params->t_0)));
      }else{
        Gnuc=1.0;
      }
    }
    return Gnuc;
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

  //---compute and store slip to determine the magnitude of an earthquake ---
  //    to this end, here the slip is computed and averaged per element
  //    in calc_seissol.f90 this value will be multiplied by the element surface
  //    and an output happened once at the end of the simulation
  void calcAverageSlip(
      std::array<real, numOfPointsPadded> &tmpSlip,
      unsigned int face
  ){
    real sum_tmpSlip = 0;
    if (m_Params->IsMagnitudeOutputOn) {
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
                         real DeltaT[CONVERGENCE_ORDER]) = 0;
};


class seissol::dr::fr_law::SolverNoFaultFL0 : public seissol::dr::fr_law::BaseFrictionSolver {

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

      //compute stresses from Qinterpolated
      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
          faultStresses.TractionGP_XY[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP];
          faultStresses.TractionGP_XZ[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP];
        }
      }
      //save stresses in imposedState
      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);
    }//End of Loop over Faces
  }//End of Function evaluate
};



class seissol::dr::fr_law::Solver_FL_33 : public seissol::dr::fr_law::BaseFrictionSolver {
protected:
  //Attributes
  real  (*nucleationStressInFaultCS)[numOfPointsPadded][6];
  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup) override {
    //first copy all Variables from the Base Lts dynRup tree
    BaseFrictionSolver::copyLtsTreeToLocal(layerData, dynRup);
    //TODO: change later to const_cast
    seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);
    nucleationStressInFaultCS =  layerData.var(ConcreteLts->nucleationStressInFaultCS); ;
    averaged_Slip             = layerData.var(ConcreteLts->averaged_Slip);
    /*
     * Add new LTS parameter specific for this
     */
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
        std::array<real, numOfPointsPadded> tmpSlip{0};
        real tn = fullUpdateTime;
        real time_inc;
        real Gnuc = 0.0;

        //compute stresses from Qinterpolated
        precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

        for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
          time_inc = DeltaT[iTimeGP];
          tn = tn + time_inc;
          Gnuc = Calc_SmoothStepIncrement(tn, time_inc)/time_inc;

          for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
            //!EQN%NucleationStressInFaultCS (1 and 2) contains the slip in FaultCS
            faultStresses.TractionGP_XY[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP] - impAndEta[ltsFace].eta_s * nucleationStressInFaultCS[ltsFace][iBndGP][0] *Gnuc;
            faultStresses.TractionGP_XZ[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP] - impAndEta[ltsFace].eta_s * nucleationStressInFaultCS[ltsFace][iBndGP][1] *Gnuc;
            slipRate1[ltsFace][iBndGP] = nucleationStressInFaultCS[ltsFace][iBndGP][0] * Gnuc;
            slipRate2[ltsFace][iBndGP] = nucleationStressInFaultCS[ltsFace][iBndGP][1] * Gnuc;
            LocSlipRate[iBndGP]  = std::sqrt( seissol::dr::aux::power(slipRate1[ltsFace][iBndGP],2) + seissol::dr::aux::power(slipRate2[ltsFace][iBndGP],2));

            //! Update slip
            slip1[ltsFace][iBndGP] += slipRate1[ltsFace][iBndGP]*time_inc;
            slip2[ltsFace][iBndGP] += slipRate2[ltsFace][iBndGP]*time_inc;
            slip[ltsFace][iBndGP] += LocSlipRate[iBndGP]*time_inc;
            tmpSlip[iBndGP] += LocSlipRate[iBndGP]*time_inc;

            tracXY[ltsFace][iBndGP] = faultStresses.TractionGP_XY[iTimeGP][iBndGP];
            tracXZ[ltsFace][iBndGP] = faultStresses.TractionGP_XY[iTimeGP][iBndGP];
          }
        }
        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        outputRuptureFront(LocSlipRate, fullUpdateTime, ltsFace);

        //output peak slip rate
        calcPeakSlipRate(LocSlipRate, ltsFace);

        //---compute and store slip to determine the magnitude of an earthquake ---
        //    to this end, here the slip is computed and averaged per element
        //    in calc_seissol.f90 this value will be multiplied by the element surface
        //    and an output happened once at the end of the simulation
        calcAverageSlip(tmpSlip, ltsFace);

        //save stresses in imposedState
        postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);
      }//End of Loop over Faces
    }
};


class seissol::dr::fr_law::SolverBluePrint : public seissol::dr::fr_law::BaseFrictionSolver {
protected:
  //Attributes
  real  (*templateAttribute)[numOfPointsPadded];

  /*
 * copies all parameters from the DynamicRupture LTS to the local attributes
 */
  void copyLtsTreeToLocal(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup) override {
    //first copy all Variables from the Base Lts dynRup tree
    BaseFrictionSolver::copyLtsTreeToLocal(layerData, dynRup);
    //TODO: change later to const_cast
    //seissol::initializers::DR_lts_template *ConcreteLts = dynamic_cast<seissol::initializers::DR_lts_template *>(dynRup);

    /*
     * Add new LTS parameter specific for this
     */
  }

  void calcSlipRate(
      FaultStresses &faultStresses,
      real LocSlipRate[seissol::tensor::resamplePar::size()],
      unsigned int iTimeGP,
      unsigned int face
  ) {
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      LocSlipRate[iBndGP] = 0;
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
      real LocSlipRate[seissol::tensor::resamplePar::size()];

      //compute stresses from Qinterpolated
      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);


      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        /*
         * add friction law calculation here:
         */
        calcSlipRate(faultStresses, LocSlipRate, iTimeGP , ltsFace);
      }
      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(LocSlipRate, fullUpdateTime, ltsFace);

      //output peak slip rate
      calcPeakSlipRate(LocSlipRate, ltsFace);

      //save stresses in imposedState
      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);

    }//End of Loop over Faces

  }//End of Function evaluate
};



#endif //SEISSOL_DR_SOLVER_BASE_H

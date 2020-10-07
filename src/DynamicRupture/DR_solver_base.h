//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_SOLVER_BASE_H
#define SEISSOL_DR_SOLVER_BASE_H

#include "DR_math.h"
#include <yaml-cpp/yaml.h>
#include <Kernels/DynamicRupture.h>


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
  dr::DrParameterT        *m_Params;
  ImpedancesAndEta*       impAndEta;
  real                    m_fullUpdateTime;
  real                    deltaT[CONVERGENCE_ORDER] = {};
  real                    (*initialStressInFaultCS)[numOfPointsPadded][6];
  real                    (*cohesion)[numOfPointsPadded];
  real                    (*mu)[numOfPointsPadded];
  real                    (*slip)[numOfPointsPadded];
  real                    (*slip1)[numOfPointsPadded];
  real                    (*slip2)[numOfPointsPadded];
  real                    (*locSlipRate)[numOfPointsPadded];
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
                                seissol::initializers::DynamicRupture *dynRup,
                                real fullUpdateTime){
    impAndEta                                     = layerData.var(dynRup->impAndEta);
    initialStressInFaultCS                        = layerData.var(dynRup->initialStressInFaultCS);
    cohesion                                      = layerData.var(dynRup->cohesion);
    mu                                            = layerData.var(dynRup->mu);
    slip                                          = layerData.var(dynRup->slip);
    slip1                                         = layerData.var(dynRup->slip1);
    slip2                                         = layerData.var(dynRup->slip2);
    locSlipRate                                   = layerData.var(dynRup->locSlipRate);
    slipRate1                                     = layerData.var(dynRup->slipRate1);
    slipRate2                                     = layerData.var(dynRup->slipRate2);
    rupture_time                                  = layerData.var(dynRup->rupture_time);
    RF                                            = layerData.var(dynRup->RF);
    peakSR                                        = layerData.var(dynRup->peakSR);
    tracXY                                        = layerData.var(dynRup->tracXY);
    tracXZ                                        = layerData.var(dynRup->tracXZ);
    imposedStatePlus                              = layerData.var(dynRup->imposedStatePlus);
    imposedStateMinus                             = layerData.var(dynRup->imposedStateMinus);
    m_fullUpdateTime                              = fullUpdateTime;
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
  virtual void precomputeStressFromQInterpolated(
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
  real Calc_SmoothStepIncrement(real current_time ,real dt){
    real Gnuc;
    real prevtime;
    if(current_time > 0.0 && current_time <= m_Params->t_0){
      Gnuc = Calc_SmoothStep(current_time);
      prevtime = current_time - dt;
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
  real Calc_SmoothStep(real current_time){
    real Gnuc;
    if (current_time <= 0){
      Gnuc=0.0;
    }else{
      if (current_time < m_Params->t_0){
        Gnuc = std::exp(seissol::dr::aux::power(current_time - m_Params->t_0, 2) / (current_time * (current_time - 2.0 * m_Params->t_0)));
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
      unsigned int ltsFace
  ){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (RF[ltsFace][iBndGP] && locSlipRate[ltsFace][iBndGP] > 0.001) {
        rupture_time[ltsFace][iBndGP] = m_fullUpdateTime;
        RF[ltsFace][iBndGP] = false;
      }
    }
  }


  void calcPeakSlipRate(
      unsigned int ltsFace){
    for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
      if (locSlipRate[ltsFace][iBndGP] > peakSR[ltsFace][iBndGP]) {
        peakSR[ltsFace][iBndGP] = locSlipRate[ltsFace][iBndGP];
      }
    }
  }

  //---compute and store slip to determine the magnitude of an earthquake ---
  //    to this end, here the slip is computed and averaged per element
  //    in calc_seissol.f90 this value will be multiplied by the element surface
  //    and an output happened once at the end of the simulation
  void calcAverageSlip(
      std::array<real, numOfPointsPadded> &tmpSlip,
      unsigned int ltsFace
  ){
    real sum_tmpSlip = 0;
    if (m_Params->IsMagnitudeOutputOn) {
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++)
        sum_tmpSlip += tmpSlip[iBndGP];
      averaged_Slip[ltsFace] = averaged_Slip[ltsFace] + sum_tmpSlip / numberOfPoints;
    }
  }

public:
  /*
   * evaluates the current friction model
   */
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                         seissol::initializers::DynamicRupture *dynRup,
                         real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                         real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                         real fullUpdateTime,
                         real timeWeights[CONVERGENCE_ORDER]) = 0;

  /*
   * compute the DeltaT from the current timePoints
   * call this function before evaluate to set the correct DeltaT
   */
  void computeDeltaT(double timePoints[CONVERGENCE_ORDER]){
    deltaT[0]= timePoints[0];
    for(int iTimeGP = 1; iTimeGP< CONVERGENCE_ORDER; iTimeGP++ ){
      deltaT[iTimeGP] = timePoints[iTimeGP]- timePoints[iTimeGP-1];
    }
    deltaT[CONVERGENCE_ORDER-1] = deltaT[CONVERGENCE_ORDER-1] + deltaT[0];  // to fill last segment of Gaussian integration
  }

  //TODO: maybe this function is used, if not delete it
  void dynamicRuptureCalc(seissol::initializers::Layer&  layerData,
                          seissol::initializers::DynamicRupture *dynRup,
                          real fullUpdateTime,
                          kernels::DynamicRupture &dynamicRuptureKernel,
                          GlobalData const* globalData){

    DRFaceInformation*                    faceInformation                                                   = layerData.var(dynRup->faceInformation);
    DRGodunovData*                        godunovData                                                       = layerData.var(dynRup->godunovData);
    real**                                timeDerivativePlus                                                = layerData.var(dynRup->timeDerivativePlus);
    real**                                timeDerivativeMinus                                               = layerData.var(dynRup->timeDerivativeMinus);
    //TODO: if computation loop is not splitted -> smaller size
    alignas(ALIGNMENT) real QInterpolatedPlus[layerData.getNumberOfCells()][CONVERGENCE_ORDER][tensor::QInterpolated::size()];
    alignas(ALIGNMENT) real QInterpolatedMinus[layerData.getNumberOfCells()][CONVERGENCE_ORDER][tensor::QInterpolated::size()];


    computeDeltaT(dynamicRuptureKernel.timePoints);

    //this copies all lts data pointers to local class attributes
    copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) //private(QInterpolatedPlus,QInterpolatedMinus)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      //initialize struct for in/outputs stresses
      FaultStresses faultStresses{};

      unsigned prefetchFace = (ltsFace < layerData.getNumberOfCells() - 1) ? ltsFace + 1 : ltsFace;
      dynamicRuptureKernel.spaceTimeInterpolation(  faceInformation[ltsFace],
                                                      globalData,
                                                      &godunovData[ltsFace],
                                                      timeDerivativePlus[ltsFace],
                                                      timeDerivativeMinus[ltsFace],
                                                      QInterpolatedPlus[ltsFace],
                                                      QInterpolatedMinus[ltsFace],
                                                      timeDerivativePlus[prefetchFace],
                                                      timeDerivativeMinus[prefetchFace] );

      //compute stresses from Qinterpolated
      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        /*
         * add friction law calculation here:
         */

      }

      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(ltsFace);

      //output peak slip rate
      calcPeakSlipRate(ltsFace);

      //save stresses in imposedState
      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, dynamicRuptureKernel.timeWeights, ltsFace);

      //*/
    } //End layerData.getNumberOfCells()-loop
  }
};  //End BaseFrictionSolver Class


class seissol::dr::fr_law::SolverNoFaultFL0 : public seissol::dr::fr_law::BaseFrictionSolver {

public:
  virtual void evaluate(seissol::initializers::Layer&  layerData,
                        seissol::initializers::DynamicRupture *dynRup,
                        real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                        real fullUpdateTime,
                        real timeWeights[CONVERGENCE_ORDER]) override {
    copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
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
                          seissol::initializers::DynamicRupture *dynRup,
                          real fullUpdateTime) override {
    //first copy all Variables from the Base Lts dynRup tree
    BaseFrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
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
                          real timeWeights[CONVERGENCE_ORDER]) override {

      copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
      for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
        //initialize struct for in/outputs stresses
        FaultStresses faultStresses{};

        //declare local variables
        std::array<real, numOfPointsPadded> tmpSlip{0};
        real tn = fullUpdateTime;
        real time_inc;
        real Gnuc = 0.0;

        //compute stresses from Qinterpolated
        precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

        for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
          time_inc = deltaT[iTimeGP];
          tn = tn + time_inc;
          Gnuc = Calc_SmoothStepIncrement(tn, time_inc)/time_inc;

          for (int iBndGP = 0; iBndGP < numberOfPoints; iBndGP++) {
            //!EQN%NucleationStressInFaultCS (1 and 2) contains the slip in FaultCS
            faultStresses.TractionGP_XY[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP] - impAndEta[ltsFace].eta_s * nucleationStressInFaultCS[ltsFace][iBndGP][0] *Gnuc;
            faultStresses.TractionGP_XZ[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP] - impAndEta[ltsFace].eta_s * nucleationStressInFaultCS[ltsFace][iBndGP][1] *Gnuc;
            slipRate1[ltsFace][iBndGP] = nucleationStressInFaultCS[ltsFace][iBndGP][0] * Gnuc;
            slipRate2[ltsFace][iBndGP] = nucleationStressInFaultCS[ltsFace][iBndGP][1] * Gnuc;
            locSlipRate[ltsFace][iBndGP]  = std::sqrt( seissol::dr::aux::power(slipRate1[ltsFace][iBndGP],2) + seissol::dr::aux::power(slipRate2[ltsFace][iBndGP],2));

            //! Update slip
            slip1[ltsFace][iBndGP] += slipRate1[ltsFace][iBndGP]*time_inc;
            slip2[ltsFace][iBndGP] += slipRate2[ltsFace][iBndGP]*time_inc;
            slip[ltsFace][iBndGP] += locSlipRate[ltsFace][iBndGP]*time_inc;
            tmpSlip[iBndGP] += locSlipRate[ltsFace][iBndGP]*time_inc;

            tracXY[ltsFace][iBndGP] = faultStresses.TractionGP_XY[iTimeGP][iBndGP];
            tracXZ[ltsFace][iBndGP] = faultStresses.TractionGP_XY[iTimeGP][iBndGP];
          }
        }
        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        outputRuptureFront(ltsFace);

        //output peak slip rate
        calcPeakSlipRate(ltsFace);

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
                          seissol::initializers::DynamicRupture *dynRup, real fullUpdateTime) override {
    //first copy all Variables from the Base Lts dynRup tree
    BaseFrictionSolver::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
    //TODO: change later to const_cast
    //seissol::initializers::DR_lts_template *ConcreteLts = dynamic_cast<seissol::initializers::DR_lts_template *>(dynRup);

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
                        real timeWeights[CONVERGENCE_ORDER]) override {
    copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      //initialize struct for in/outputs stresses
      FaultStresses faultStresses{};

      //compute stresses from Qinterpolated
      precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);


      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
        /*
         * add friction law calculation here:
         */
      }
      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      outputRuptureFront(ltsFace);

      //output peak slip rate
      calcPeakSlipRate(ltsFace);

      //save stresses in imposedState
      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);

    }//End of Loop over Faces

  }//End of Function evaluate
};



#endif //SEISSOL_DR_SOLVER_BASE_H

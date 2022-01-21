#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Parameters.h"
#include "FrictionSolver.h"

namespace seissol::dr::friction_law {
/**
 * Base class, has implementations of methods that are used by each friction law
 * Actual friction law is plugged in via CRTP.
 */
template <typename Derived>
class BaseFrictionLaw : public FrictionSolver {
  public:
  BaseFrictionLaw(dr::DRParameters& drParameters) : drParameters(drParameters){};

  dr::DRParameters& drParameters;
  ImpedancesAndEta* impAndEta;
  real mFullUpdateTime;
  // CS = coordinate system
  real (*initialStressInFaultCS)[misc::numPaddedPoints][6];
  real (*cohesion)[misc::numPaddedPoints];
  real (*mu)[misc::numPaddedPoints];
  real (*accumulatedSlipMagnitude)[misc::numPaddedPoints];
  real (*slip1)[misc::numPaddedPoints];
  real (*slip2)[misc::numPaddedPoints];
  real (*slipRateMagnitude)[misc::numPaddedPoints];
  real (*slipRate1)[misc::numPaddedPoints];
  real (*slipRate2)[misc::numPaddedPoints];
  real (*ruptureTime)[misc::numPaddedPoints];
  bool (*ruptureTimePending)[misc::numPaddedPoints];
  real (*peakSlipRate)[misc::numPaddedPoints];
  real (*tractionXY)[misc::numPaddedPoints];
  real (*tractionXZ)[misc::numPaddedPoints];
  real (*imposedStatePlus)[tensor::QInterpolated::size()];
  real (*imposedStateMinus)[tensor::QInterpolated::size()];

  // be careful only for some FLs initialized:
  real* averagedSlip;
  real (*dynStressTime)[misc::numPaddedPoints];
  bool (*dynStressTimePending)[misc::numPaddedPoints];

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    impAndEta = layerData.var(dynRup->impAndEta);
    initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS);
    mu = layerData.var(dynRup->mu);
    accumulatedSlipMagnitude = layerData.var(dynRup->accumulatedSlipMagnitude);
    slip1 = layerData.var(dynRup->slip1);
    slip2 = layerData.var(dynRup->slip2);
    slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude);
    slipRate1 = layerData.var(dynRup->slipRate1);
    slipRate2 = layerData.var(dynRup->slipRate2);
    ruptureTime = layerData.var(dynRup->ruptureTime);
    ruptureTimePending = layerData.var(dynRup->ruptureTimePending);
    peakSlipRate = layerData.var(dynRup->peakSlipRate);
    tractionXY = layerData.var(dynRup->tractionXY);
    tractionXZ = layerData.var(dynRup->tractionXZ);
    imposedStatePlus = layerData.var(dynRup->imposedStatePlus);
    imposedStateMinus = layerData.var(dynRup->imposedStateMinus);
    mFullUpdateTime = fullUpdateTime;
    averagedSlip = layerData.var(dynRup->averagedSlip);
    dynStressTime = layerData.var(dynRup->dynStressTime);
    dynStressTimePending = layerData.var(dynRup->dynStressTimePending);
    static_cast<Derived*>(this)->copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
  }

  /**
   * Calculate godunov state from jump of plus and minus side
   * using equations (A2) from Pelites et al. 2014
   * Definiton of eta and impedance Z are found in dissertation of Carsten Uphoff
   * input:
   * QInterpolatedPlus, QInterpolatedMinus
   * output:
   * NorStressGP, XYStressGP, XZStressGP
   */
  void precomputeStressFromQInterpolated(
      FaultStresses& faultStresses,
      real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      unsigned int ltsFace) {

    static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                  "Different number of quadrature points?");

    // this initialization of the kernel could be moved to the initializer,
    // since all inputs outside the j-loop are time independent
    // set inputParam could be extendent for this
    // the kernel then could be a class attribute (but be careful of race conditions since this is
    // computed in parallel!!)
    dynamicRupture::kernel::StressFromQInterpolated stressFromQInterpolatedKrnl;
    stressFromQInterpolatedKrnl.eta_p = impAndEta[ltsFace].etaP;
    stressFromQInterpolatedKrnl.eta_s = impAndEta[ltsFace].etaS;
    stressFromQInterpolatedKrnl.inv_Zp = impAndEta[ltsFace].invZp;
    stressFromQInterpolatedKrnl.inv_Zs = impAndEta[ltsFace].invZs;
    stressFromQInterpolatedKrnl.inv_Zp_neig = impAndEta[ltsFace].invZpNeig;
    stressFromQInterpolatedKrnl.inv_Zs_neig = impAndEta[ltsFace].invZsNeig;
    stressFromQInterpolatedKrnl.select0 = init::select0::Values;
    stressFromQInterpolatedKrnl.select3 = init::select3::Values;
    stressFromQInterpolatedKrnl.select5 = init::select5::Values;
    stressFromQInterpolatedKrnl.select6 = init::select6::Values;
    stressFromQInterpolatedKrnl.select7 = init::select7::Values;
    stressFromQInterpolatedKrnl.select8 = init::select8::Values;

    for (unsigned j = 0; j < CONVERGENCE_ORDER; j++) {
      stressFromQInterpolatedKrnl.QInterpolatedMinus = qInterpolatedMinus[j];
      stressFromQInterpolatedKrnl.QInterpolatedPlus = qInterpolatedPlus[j];
      stressFromQInterpolatedKrnl.NorStressGP = faultStresses.normalStress[j];
      stressFromQInterpolatedKrnl.XYStressGP = faultStresses.xyStress[j];
      stressFromQInterpolatedKrnl.XZStressGP = faultStresses.xzStress[j];
      // Carsten Uphoff Thesis: EQ.: 4.53
      stressFromQInterpolatedKrnl.execute();
    }
  }

  /**
   * Integrate over all Time points with the time weights and calculate the traction vor each side
   * according to Carsten Uphoff Thesis: EQ.: 4.60
   * IN: NormalStressGP, XYTractionResultGP, * XZTractionResultGP
   * OUT: imposedStatePlus, imposedStateMinus
   */
  void postcomputeImposedStateFromNewStress(
      real qInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real qInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      const FaultStresses& faultStresses,
      double timeWeights[CONVERGENCE_ORDER],
      unsigned int ltsFace) {
    // this initialization of the kernel could be moved to the initializer
    // set inputParam could be extendent for this (or create own function)
    // the kernel then could be a class attribute and following values are only set once
    //(but be careful of race conditions since this is computed in parallel for each face!!)
    dynamicRupture::kernel::ImposedStateFromNewStress imposedStateFromNewStressKrnl;
    imposedStateFromNewStressKrnl.select0 = init::select0::Values;
    imposedStateFromNewStressKrnl.select3 = init::select3::Values;
    imposedStateFromNewStressKrnl.select5 = init::select5::Values;
    imposedStateFromNewStressKrnl.select6 = init::select6::Values;
    imposedStateFromNewStressKrnl.select7 = init::select7::Values;
    imposedStateFromNewStressKrnl.select8 = init::select8::Values;
    imposedStateFromNewStressKrnl.inv_Zs = impAndEta[ltsFace].invZs;
    imposedStateFromNewStressKrnl.inv_Zs_neig = impAndEta[ltsFace].invZsNeig;
    imposedStateFromNewStressKrnl.inv_Zp = impAndEta[ltsFace].invZp;
    imposedStateFromNewStressKrnl.inv_Zp_neig = impAndEta[ltsFace].invZpNeig;

    // set imposed state to zero
    for (unsigned int i = 0; i < tensor::QInterpolated::size(); i++) {
      imposedStatePlus[ltsFace][i] = 0;
      imposedStateMinus[ltsFace][i] = 0;
    }
    imposedStateFromNewStressKrnl.imposedStatePlus = imposedStatePlus[ltsFace];
    imposedStateFromNewStressKrnl.imposedStateMinus = imposedStateMinus[ltsFace];

    for (unsigned j = 0; j < CONVERGENCE_ORDER; j++) {
      imposedStateFromNewStressKrnl.NorStressGP = faultStresses.normalStress[j];
      imposedStateFromNewStressKrnl.TractionGP_XY = faultStresses.xyTractionResult[j];
      imposedStateFromNewStressKrnl.TractionGP_XZ = faultStresses.xzTractionResult[j];
      imposedStateFromNewStressKrnl.timeWeights = timeWeights[j];
      imposedStateFromNewStressKrnl.QInterpolatedMinus = qInterpolatedMinus[j];
      imposedStateFromNewStressKrnl.QInterpolatedPlus = qInterpolatedPlus[j];
      // Carsten Uphoff Thesis: EQ.: 4.60
      imposedStateFromNewStressKrnl.execute();
    }
  }

  /**
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
   */
  real calcSmoothStepIncrement(real currentTime, real dt) {
    real gNuc = calcSmoothStep(currentTime);
    real prevTime = currentTime - dt;
    gNuc = gNuc - calcSmoothStep(prevTime);
    return gNuc;
  }

  /**
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
   */
  real calcSmoothStep(real currentTime) {
    if (currentTime <= 0) {
      return 0.0;
    } else if (currentTime < drParameters.t0) {
      return std::exp(misc::power<2>(currentTime - drParameters.t0) /
                      (currentTime * (currentTime - 2.0 * drParameters.t0)));
    } else {
      return 1.0;
    }
  }

  /**
   * output rupture front, saves update time of the rupture front
   * rupture front is the first registered change in slip rates that exceeds 0.001
   */
  void saveRuptureFrontOutput(unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      constexpr real ruptureFrontThreshold = 0.001;
      if (ruptureTimePending[ltsFace][pointIndex] &&
          slipRateMagnitude[ltsFace][pointIndex] > ruptureFrontThreshold) {
        ruptureTime[ltsFace][pointIndex] = mFullUpdateTime;
        ruptureTimePending[ltsFace][pointIndex] = false;
      }
    }
  }

  /**
   * Save the maximal computed slip rate magnitude in peakSlipRate
   */
  void savePeakSlipRateOutput(unsigned int ltsFace) {
    for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
      peakSlipRate[ltsFace][pointIndex] =
          std::max(peakSlipRate[ltsFace][pointIndex], slipRateMagnitude[ltsFace][pointIndex]);
    }
  }

  /**
   * Compute and store element-averaged slip to determine the magnitude of an earthquake.
   * In calc_seissol.f90 this value will be multiplied by the element surface
   * and the seismic moment is outputted once at the end of the simulation.
   * @param tmpSlip
   * @param ltsFace
   */
  void saveAverageSlipOutput(std::array<real, misc::numPaddedPoints>& tmpSlip,
                             unsigned int ltsFace) {
    real sumOfTmpSlip = 0;
    if (drParameters.isMagnitudeOutputOn) {
      for (unsigned pointIndex = 0; pointIndex < misc::numberOfBoundaryGaussPoints; pointIndex++) {
        sumOfTmpSlip += tmpSlip[pointIndex];
      }
      averagedSlip[ltsFace] += sumOfTmpSlip / misc::numberOfBoundaryGaussPoints;
    }
  }

  /**
   * evaluates the current friction model
   */
  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture* dynRup,
                real (*qInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real (*qInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real fullUpdateTime,
                double timeWeights[CONVERGENCE_ORDER]) {
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // loop over all dynamic rupture faces, in this LTS layer
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      // initialize struct for in/outputs stresses
      FaultStresses faultStresses = {};

      this->precomputeStressFromQInterpolated(
          faultStresses, qInterpolatedPlus[ltsFace], qInterpolatedMinus[ltsFace], ltsFace);

      // define some temporary variables
      std::array<real, misc::numPaddedPoints> stateVariableBuffer{0};
      std::array<real, misc::numPaddedPoints> strengthBuffer{0};

      static_cast<Derived*>(this)->preHook(stateVariableBuffer, ltsFace);

      // loop over sub time steps (i.e. quadrature points in time)
      for (unsigned timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
        static_cast<Derived*>(this)->updateFrictionAndSlip(
            faultStresses, stateVariableBuffer, strengthBuffer, ltsFace, timeIndex);
      }

      static_cast<Derived*>(this)->postHook(stateVariableBuffer, ltsFace);

      // output rupture front
      this->saveRuptureFrontOutput(ltsFace);

      // output time when shear stress is equal to the dynamic stress after rupture arrived
      static_cast<Derived*>(this)->saveDynamicStressOutput(ltsFace);

      // output peak slip rate
      this->savePeakSlipRateOutput(ltsFace);

      // output average slip
      // TODO: What about outputSlip
      // this->saveAverageSlipOutput(outputSlip, ltsFace);

      // compute output
      this->postcomputeImposedStateFromNewStress(qInterpolatedPlus[ltsFace],
                                                 qInterpolatedMinus[ltsFace],
                                                 faultStresses,
                                                 timeWeights,
                                                 ltsFace);
    }
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H

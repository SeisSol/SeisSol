#ifndef SEISSOL_BASEFRICTIONLAW_H
#define SEISSOL_BASEFRICTIONLAW_H

#include <yaml-cpp/yaml.h>

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
  real m_fullUpdateTime;
  // CS = coordinate system
  real (*initialStressInFaultCS)[numPaddedPoints][6];
  real (*cohesion)[numPaddedPoints];
  real (*mu)[numPaddedPoints];
  real (*slip)[numPaddedPoints];
  real (*slipStrike)[numPaddedPoints];
  real (*slipDip)[numPaddedPoints];
  real (*slipRateMagnitude)[numPaddedPoints];
  real (*slipRateStrike)[numPaddedPoints];
  real (*slipRateDip)[numPaddedPoints];
  real (*ruptureTime)[numPaddedPoints];
  bool (*ruptureFront)[numPaddedPoints];
  real (*peakSlipRate)[numPaddedPoints];
  real (*tractionXY)[numPaddedPoints];
  real (*tractionXZ)[numPaddedPoints];
  real (*imposedStatePlus)[tensor::QInterpolated::size()];
  real (*imposedStateMinus)[tensor::QInterpolated::size()];

  // be careful only for some FLs initialized:
  real* averagedSlip;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) {
    impAndEta = layerData.var(dynRup->impAndEta);
    initialStressInFaultCS = layerData.var(dynRup->initialStressInFaultCS);
    mu = layerData.var(dynRup->mu);
    slip = layerData.var(dynRup->slip);
    slipStrike = layerData.var(dynRup->slipStrike);
    slipDip = layerData.var(dynRup->slipDip);
    slipRateMagnitude = layerData.var(dynRup->slipRateMagnitude);
    slipRateStrike = layerData.var(dynRup->slipRateStrike);
    slipRateDip = layerData.var(dynRup->slipRateDip);
    ruptureTime = layerData.var(dynRup->ruptureTime);
    ruptureFront = layerData.var(dynRup->ruptureFront);
    peakSlipRate = layerData.var(dynRup->peakSlipRate);
    tractionXY = layerData.var(dynRup->tractionXY);
    tractionXZ = layerData.var(dynRup->tractionXZ);
    imposedStatePlus = layerData.var(dynRup->imposedStatePlus);
    imposedStateMinus = layerData.var(dynRup->imposedStateMinus);
    m_fullUpdateTime = fullUpdateTime;
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
      real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      unsigned int ltsFace) {

    static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0],
                  "Different number of quadrature points?");

    // this initialization of the kernel could be moved to the initializer,
    // since all inputs outside the j-loop are time independent
    // set inputParam could be extendent for this
    // the kernel then could be a class attribute (but be careful of race conditions since this is
    // computed in parallel!!)
    dynamicRupture::kernel::StressFromQInterpolated StressFromQInterpolatedKrnl;
    StressFromQInterpolatedKrnl.eta_p = impAndEta[ltsFace].eta_p;
    StressFromQInterpolatedKrnl.eta_s = impAndEta[ltsFace].eta_s;
    StressFromQInterpolatedKrnl.inv_Zp = impAndEta[ltsFace].inv_Zp;
    StressFromQInterpolatedKrnl.inv_Zs = impAndEta[ltsFace].inv_Zs;
    StressFromQInterpolatedKrnl.inv_Zp_neig = impAndEta[ltsFace].inv_Zp_neig;
    StressFromQInterpolatedKrnl.inv_Zs_neig = impAndEta[ltsFace].inv_Zs_neig;
    StressFromQInterpolatedKrnl.select0 = init::select0::Values;
    StressFromQInterpolatedKrnl.select3 = init::select3::Values;
    StressFromQInterpolatedKrnl.select5 = init::select5::Values;
    StressFromQInterpolatedKrnl.select6 = init::select6::Values;
    StressFromQInterpolatedKrnl.select7 = init::select7::Values;
    StressFromQInterpolatedKrnl.select8 = init::select8::Values;

    for (int j = 0; j < CONVERGENCE_ORDER; j++) {
      StressFromQInterpolatedKrnl.QInterpolatedMinus = QInterpolatedMinus[j];
      StressFromQInterpolatedKrnl.QInterpolatedPlus = QInterpolatedPlus[j];
      StressFromQInterpolatedKrnl.NorStressGP = faultStresses.NormalStressGP[j];
      StressFromQInterpolatedKrnl.XYStressGP = faultStresses.XYStressGP[j];
      StressFromQInterpolatedKrnl.XZStressGP = faultStresses.XZStressGP[j];
      // Carsten Uphoff Thesis: EQ.: 4.53
      StressFromQInterpolatedKrnl.execute();
    }
  }

  /**
   * Integrate over all Time points with the time weights and calculate the traction vor each side
   * according to Carsten Uphoff Thesis: EQ.: 4.60
   * IN: NormalStressGP, XYTractionResultGP, * XZTractionResultGP
   * OUT: imposedStatePlus, imposedStateMinus
   */
  void postcomputeImposedStateFromNewStress(
      real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
      const FaultStresses& faultStresses,
      double timeWeights[CONVERGENCE_ORDER],
      unsigned int ltsFace) {
    // this initialization of the kernel could be moved to the initializer
    // set inputParam could be extendent for this (or create own function)
    // the kernel then could be a class attribute and following values are only set once
    //(but be careful of race conditions since this is computed in parallel for each face!!)
    dynamicRupture::kernel::ImposedStateFromNewStress ImposedStateFromNewStressKrnl;
    ImposedStateFromNewStressKrnl.select0 = init::select0::Values;
    ImposedStateFromNewStressKrnl.select3 = init::select3::Values;
    ImposedStateFromNewStressKrnl.select5 = init::select5::Values;
    ImposedStateFromNewStressKrnl.select6 = init::select6::Values;
    ImposedStateFromNewStressKrnl.select7 = init::select7::Values;
    ImposedStateFromNewStressKrnl.select8 = init::select8::Values;
    ImposedStateFromNewStressKrnl.inv_Zs = impAndEta[ltsFace].inv_Zs;
    ImposedStateFromNewStressKrnl.inv_Zs_neig = impAndEta[ltsFace].inv_Zs_neig;
    ImposedStateFromNewStressKrnl.inv_Zp = impAndEta[ltsFace].inv_Zp;
    ImposedStateFromNewStressKrnl.inv_Zp_neig = impAndEta[ltsFace].inv_Zp_neig;

    // set imposed state to zero
    for (unsigned int i = 0; i < tensor::QInterpolated::size(); i++) {
      imposedStatePlus[ltsFace][i] = 0;
      imposedStateMinus[ltsFace][i] = 0;
    }
    ImposedStateFromNewStressKrnl.imposedStatePlus = imposedStatePlus[ltsFace];
    ImposedStateFromNewStressKrnl.imposedStateMinus = imposedStateMinus[ltsFace];

    for (int j = 0; j < CONVERGENCE_ORDER; j++) {
      ImposedStateFromNewStressKrnl.NorStressGP = faultStresses.NormalStressGP[j];
      ImposedStateFromNewStressKrnl.TractionGP_XY = faultStresses.XYTractionResultGP[j];
      ImposedStateFromNewStressKrnl.TractionGP_XZ = faultStresses.XZTractionResultGP[j];
      ImposedStateFromNewStressKrnl.timeWeights = timeWeights[j];
      ImposedStateFromNewStressKrnl.QInterpolatedMinus = QInterpolatedMinus[j];
      ImposedStateFromNewStressKrnl.QInterpolatedPlus = QInterpolatedPlus[j];
      // Carsten Uphoff Thesis: EQ.: 4.60
      ImposedStateFromNewStressKrnl.execute();
    }
  }

  /**
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
   */
  real calcSmoothStepIncrement(real currentTime, real dt) {
    if (currentTime > 0.0 && currentTime <= drParameters.t0) {
      real gNuc = calcSmoothStep(currentTime);
      real prevTime = currentTime - dt;
      if (prevTime > 0.0) {
        gNuc = gNuc - calcSmoothStep(prevTime);
      }
      return gNuc;
    } else {
      return 0.0;
    }
  }

  /**
   * For reference, see: https://strike.scec.org/cvws/download/SCEC_validation_slip_law.pdf
   */
  real calcSmoothStep(real currentTime) {
    if (currentTime <= 0) {
      return 0.0;
    } else if (currentTime < drParameters.t0) {
      return std::exp(std::pow(currentTime - drParameters.t0, 2) /
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
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      constexpr real ruptureFrontThreshold = 0.001;
      if (ruptureFront[ltsFace][pointIndex] &&
          slipRateMagnitude[ltsFace][pointIndex] > ruptureFrontThreshold) {
        ruptureTime[ltsFace][pointIndex] = m_fullUpdateTime;
        ruptureFront[ltsFace][pointIndex] = false;
      }
    }
  }

  /**
   * save the maximal computed slip rate magnitude in RpeakSlipRate
   */
  void savePeakSlipRateOutput(unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      if (slipRateMagnitude[ltsFace][pointIndex] > peakSlipRate[ltsFace][pointIndex]) {
        peakSlipRate[ltsFace][pointIndex] = slipRateMagnitude[ltsFace][pointIndex];
      }
    }
  }

  /**
   * compute and store slip to determine the magnitude of an earthquake
   * to this end, here the slip is computed and averaged per element
   * in calc_seissol.f90 this value will be multiplied by the element surface
   * and an output happened once at the end of the simulation
   * @param tmpSlip
   * @param ltsFace
   */
  void saveAverageSlipOutput(std::array<real, numPaddedPoints>& tmpSlip, unsigned int ltsFace) {
    real sumOfTmpSlip = 0;
    if (drParameters.isMagnitudeOutputOn) {
      for (int pointIndex = 0; pointIndex < numberOfPoints; pointIndex++) {
        sumOfTmpSlip += tmpSlip[pointIndex];
      }
      averagedSlip[ltsFace] = averagedSlip[ltsFace] + sumOfTmpSlip / numberOfPoints;
    }
  }

  /**
   *  compute the slip rate and the traction from the fault strength and fault stresses
   *  also updates the directional slipStrike and slipDip
   */
  void calcSlipRateAndTraction(FaultStresses& faultStresses,
                               std::array<real, numPaddedPoints>& strength,
                               unsigned int timeIndex,
                               unsigned int ltsFace) {
    for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
      // calculate absolute value of stress in Y and Z direction
      real totalStressXY = initialStressInFaultCS[ltsFace][pointIndex][3] +
                           faultStresses.XYStressGP[timeIndex][pointIndex];
      real totalStressXZ = initialStressInFaultCS[ltsFace][pointIndex][5] +
                           faultStresses.XZStressGP[timeIndex][pointIndex];
      real absoluteShearStress = std::sqrt(std::pow(totalStressXY, 2) + std::pow(totalStressXZ, 2));

      // calculate slip rates
      slipRateMagnitude[ltsFace][pointIndex] =
          std::max(static_cast<real>(0.0),
                   (absoluteShearStress - strength[pointIndex]) * impAndEta[ltsFace].inv_eta_s);

      slipRateStrike[ltsFace][pointIndex] =
          slipRateMagnitude[ltsFace][pointIndex] * totalStressXY / absoluteShearStress;
      slipRateDip[ltsFace][pointIndex] =
          slipRateMagnitude[ltsFace][pointIndex] * totalStressXZ / absoluteShearStress;

      // calculate traction
      faultStresses.XYTractionResultGP[timeIndex][pointIndex] =
          faultStresses.XYStressGP[timeIndex][pointIndex] -
          impAndEta[ltsFace].eta_s * slipRateStrike[ltsFace][pointIndex];
      faultStresses.XZTractionResultGP[timeIndex][pointIndex] =
          faultStresses.XZStressGP[timeIndex][pointIndex] -
          impAndEta[ltsFace].eta_s * slipRateDip[ltsFace][pointIndex];
      tractionXY[ltsFace][pointIndex] = faultStresses.XYTractionResultGP[timeIndex][pointIndex];
      tractionXZ[ltsFace][pointIndex] = faultStresses.XYTractionResultGP[timeIndex][pointIndex];

      // update directional slip
      slipStrike[ltsFace][pointIndex] += slipRateStrike[ltsFace][pointIndex] * deltaT[timeIndex];
      slipDip[ltsFace][pointIndex] += slipRateDip[ltsFace][pointIndex] * deltaT[timeIndex];
    }
  }

  /**
   * evaluates the current friction model
   * TODO: move implementation here and define useful hooks
   */
  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture* dynRup,
                real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
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
          faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      // define some temporary variables
      std::array<real, numPaddedPoints> stateVariableBuffer{0};
      std::array<real, numPaddedPoints> strengthBuffer{0};

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
      // Todo: Only for linear slip weakening
      // this->saveDynamicStressOutput(ltsFace);

      // output peak slip rate
      this->savePeakSlipRateOutput(ltsFace);

      // output average slip
      // TODO: What about outputSlip
      // this->saveAverageSlipOutput(outputSlip, ltsFace);

      // compute output
      this->postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                                 QInterpolatedMinus[ltsFace],
                                                 faultStresses,
                                                 timeWeights,
                                                 ltsFace);
    }
  }
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_BASEFRICTIONLAW_H

#ifndef SEISSOL_FRICTIONLAW_BLUEPRINT_H
#define SEISSOL_FRICTIONLAW_BLUEPRINT_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
class BluePrint;
};
/*
 * examplary implementation of a new friction law structure. Can be used as template to create new
 * friction laws.
 */
class seissol::dr::friction_law::BluePrint : public seissol::dr::friction_law::BaseFrictionLaw {
  protected:
  // Attributes
  real (*templateAttribute)[numOfPointsPadded];

  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) override {
    // first copy all Variables from the Base Lts dynRup tree
    BaseFrictionLaw::copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);

    // seissol::initializers::DR_lts_template *ConcreteLts =
    // dynamic_cast<seissol::initializers::DR_lts_template *>(dynRup);

    /*
     * Add new LTS parameter specific for this
     */
  }

  public:
  virtual void
      evaluate(seissol::initializers::Layer& layerData,
               seissol::initializers::DynamicRupture* dynRup,
               real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real fullUpdateTime,
               real timeWeights[CONVERGENCE_ORDER]) override {
    copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
      // initialize struct for in/outputs stresses
      FaultStresses faultStresses = {};

      // compute stresses from Qinterpolated
      precomputeStressFromQInterpolated(
          faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

      for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) { // loop over time steps
        /*
         * add friction law calculation here:
         * computed slip rates, traction, friction coefficients and state variables
         */
      }
      // output rupture front
      // outside of iTimeGP loop in order to safe an 'if' in a loop
      // this way, no subtimestep resolution possible
      saveRuptureFrontOutput(ltsFace);

      // output peak slip rate
      savePeakSlipRateOutput(ltsFace);

      // save stresses in imposedState
      postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                           QInterpolatedMinus[ltsFace],
                                           faultStresses,
                                           timeWeights,
                                           ltsFace);

    } // End of Loop over Faces

  } // End of Function evaluate
};

#endif // SEISSOL_FRICTIONLAW_BLUEPRINT_H

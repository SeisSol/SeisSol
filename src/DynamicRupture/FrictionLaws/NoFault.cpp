#include "NoFault.h"

  void seissol::dr::friction_law::NoFault::evaluate(seissol::initializers::Layer &layerData,
                                 seissol::initializers::DynamicRupture *dynRup,
                                 real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                                 real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                                 real fullUpdateTime,
                                 real timeWeights[CONVERGENCE_ORDER]) {
  copyLtsTreeToLocal(layerData, dynRup, fullUpdateTime);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned ltsFace = 0; ltsFace < layerData.getNumberOfCells(); ++ltsFace) {
    //initialize struct for in/outputs stresses
    FaultStresses faultStresses = {};

    //compute stresses from Qinterpolated
    precomputeStressFromQInterpolated(faultStresses, QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], ltsFace);

    for (int iTimeGP = 0; iTimeGP < CONVERGENCE_ORDER; iTimeGP++) {  //loop over time steps
      for (int iBndGP = 0; iBndGP < numOfPointsPadded; iBndGP++) {
        faultStresses.XYTractionResultGP[iTimeGP][iBndGP] = faultStresses.XYStressGP[iTimeGP][iBndGP];
        faultStresses.XZTractionResultGP[iTimeGP][iBndGP] = faultStresses.XZStressGP[iTimeGP][iBndGP];
      }
    }
    //save stresses in imposedState
    postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace], QInterpolatedMinus[ltsFace], faultStresses, timeWeights, ltsFace);
  }//End of Loop over Faces
}//End of Function evaluate



#include "NoFault.h"

namespace seissol::dr::friction_law {
void NoFault::evaluate(seissol::initializers::Layer& layerData,
                       seissol::initializers::DynamicRupture* dynRup,
                       real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                       real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                       real fullUpdateTime,
                       double timeWeights[CONVERGENCE_ORDER]) {
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

    for (int timeIndex = 0; timeIndex < CONVERGENCE_ORDER; timeIndex++) { // loop over time steps
      for (int pointIndex = 0; pointIndex < numPaddedPoints; pointIndex++) {
        faultStresses.XYTractionResultGP[timeIndex][pointIndex] =
            faultStresses.XYStressGP[timeIndex][pointIndex];
        faultStresses.XZTractionResultGP[timeIndex][pointIndex] =
            faultStresses.XZStressGP[timeIndex][pointIndex];
      }
    }
    // save stresses in imposedState
    postcomputeImposedStateFromNewStress(QInterpolatedPlus[ltsFace],
                                         QInterpolatedMinus[ltsFace],
                                         faultStresses,
                                         timeWeights,
                                         ltsFace);
  } // End of Loop over Faces
} // End of Function evaluate
} // namespace seissol::dr::friction_law
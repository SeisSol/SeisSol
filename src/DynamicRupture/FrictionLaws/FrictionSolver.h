#ifndef SEISSOL_FRICTIONSOLVER_H
#define SEISSOL_FRICTIONSOLVER_H

#include "DynamicRupture/Misc.h"
#include "Initializer/DynamicRupture.h"
#include "Kernels/DynamicRupture.h"

namespace seissol::dr::friction_law {
/**
 * Struct that contains all input stresses and output stresses
 * IN: NormalStress, XYStress, XZStress (Godunov stresses computed by
 * precomputeStressFromQInterpolated) OUT: XYTractionResult, XZTractionResult
 * and NormalStress (used to compute resulting +/- sided stress results by
 * postcomputeImposedStateFromNewStress)
 */
struct FaultStresses {
  real xyTractionResult[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real xzTractionResult[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real normalStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real xyStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
  real xzStress[CONVERGENCE_ORDER][misc::numPaddedPoints] = {{}};
};

/**
 * Abstract Base for friction solver class with the public interface
 * Only needed to be able to store a shared_ptr<FrictionSolver> in MemoryManager and TimeCluster.
 * BaseFrictionLaw has a template argument for CRTP, hence, we can't store a pointer to any
 * BaseFrictionLaw.
 */
class FrictionSolver {
  public:
  virtual ~FrictionSolver(){};

  virtual void
      evaluate(seissol::initializers::Layer& layerData,
               seissol::initializers::DynamicRupture* dynRup,
               real (*qInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real (*qInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real fullUpdateTime,
               double timeWeights[CONVERGENCE_ORDER]) = 0;

  /**
   * compute the DeltaT from the current timePoints
   * call this function before evaluate to set the correct DeltaT
   */
  void computeDeltaT(double timePoints[CONVERGENCE_ORDER]) {
    deltaT[0] = timePoints[0];
    for (unsigned timeIndex = 1; timeIndex < CONVERGENCE_ORDER; timeIndex++) {
      deltaT[timeIndex] = timePoints[timeIndex] - timePoints[timeIndex - 1];
    }
    // to fill last segment of Gaussian integration
    deltaT[CONVERGENCE_ORDER - 1] = deltaT[CONVERGENCE_ORDER - 1] + deltaT[0];
  }

  real deltaT[CONVERGENCE_ORDER] = {};
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_FRICTIONSOLVER_H

#ifndef SEISSOL_VELOCITYWEAKENING_H
#define SEISSOL_VELOCITYWEAKENING_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * This class was not tested and compared to the Fortran FL4. Since FL4 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class VelocityWeakening : public BaseFrictionLaw {
  protected:
  // Attributes
  real (*stateVar)[numPaddedPoints]{};
  real (*sl0)[numPaddedPoints]{};
  real (*a)[numPaddedPoints]{};

  /*
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime) override;

  public:
  virtual void
      evaluate(seissol::initializers::Layer& layerData,
               seissol::initializers::DynamicRupture* dynRup,
               real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real fullUpdateTime,
               double timeWeights[CONVERGENCE_ORDER]) override;
};
} // namespace seissol::dr::friction_law

#endif // SEISSOL_VELOCITYWEAKENING_H

#ifndef SEISSOL_IMPOSEDSLIPRATES_H
#define SEISSOL_IMPOSEDSLIPRATES_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
class ImposedSlipRates;
}
/*
 * Slip rates are set fixed values (defined by nucleationStressInFaultCS)
 */
class seissol::dr::friction_law::ImposedSlipRates : public BaseFrictionLaw {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  protected:
  // Attributes
  real (*nucleationStressInFaultCS)[numPaddedPoints][6];

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

#endif // SEISSOL_IMPOSEDSLIPRATES_H

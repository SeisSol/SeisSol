#ifndef SEISSOL_IMPOSEDSLIPRATES_H
#define SEISSOL_IMPOSEDSLIPRATES_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * Slip rates are set fixed values
 */
class ImposedSlipRates : public BaseFrictionLaw {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  protected:
  // CS = coordinate system
  real (*nucleationStressInFaultCS)[numPaddedPoints][6];

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
#endif // SEISSOL_IMPOSEDSLIPRATES_H

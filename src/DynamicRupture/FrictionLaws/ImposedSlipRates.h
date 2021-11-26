#ifndef SEISSOL_IMPOSEDSLIPRATES_H
#define SEISSOL_IMPOSEDSLIPRATES_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
/**
 * Slip rates are set fixed values
 */
class ImposedSlipRates : public BaseFrictionLaw<ImposedSlipRates> {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  // CS = coordinate system
  real (*nucleationStressInFaultCS)[numPaddedPoints][6];

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             std::array<real, numPaddedPoints>& stateVariableBuffer,
                             std::array<real, numPaddedPoints>& strengthBuffer,
                             unsigned& ltsFace,
                             unsigned& timeIndex);

  void preHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace);
  void postHook(std::array<real, numPaddedPoints>& stateVariableBuffer, unsigned ltsFace);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_IMPOSEDSLIPRATES_H

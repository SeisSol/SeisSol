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

  real (*strikeSlip)[misc::numPaddedPoints];
  real (*dipSlip)[misc::numPaddedPoints];
  real (*onsetTime)[misc::numPaddedPoints];
  real (*tauS)[misc::numPaddedPoints];
  real (*tauR)[misc::numPaddedPoints];

  void copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                          seissol::initializers::DynamicRupture* dynRup,
                          real fullUpdateTime);

  void updateFrictionAndSlip(FaultStresses& faultStresses,
                             TractionResults& tractionResults,
                             std::array<real, misc::numPaddedPoints>& stateVariableBuffer,
                             std::array<real, misc::numPaddedPoints>& strengthBuffer,
                             unsigned& ltsFace,
                             unsigned& timeIndex);

  void preHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned ltsFace);
  void postHook(std::array<real, misc::numPaddedPoints>& stateVariableBuffer, unsigned ltsFace);
  void saveDynamicStressOutput(unsigned int ltsFace);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_IMPOSEDSLIPRATES_H

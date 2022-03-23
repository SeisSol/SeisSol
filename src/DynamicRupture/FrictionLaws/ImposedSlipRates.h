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

  private:
  /**
   * Implementation of the regularized Yoffe function * defined in Appendix of Tinti et al. (2005)
   */
  real regularizedYoffe(real time, real tauS, real tauR);
  /**
   * c1 to c6 are analytical functions * used for building the regularized Yoffe function
   */
  real c1(real time, real tauS, real tauR);
  real c2(real time, real tauS, real tauR);
  real c3(real time, real tauS, real tauR);
  real c4(real time, real tauS, real tauR);
  real c5(real time, real tauS, real tauR);
  real c6(real time, real tauS, real tauR);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_IMPOSEDSLIPRATES_H

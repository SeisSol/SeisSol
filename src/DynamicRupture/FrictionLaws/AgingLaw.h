#ifndef SEISSOL_AGINGLAW_H
#define SEISSOL_AGINGLAW_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {
class AgingLaw;
}

/*
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class seissol::dr::friction_law::AgingLaw : public BaseFrictionLaw {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;

  protected:
  virtual real calcStateVariableHook(real SV0, real tmp, real time_inc, real RS_sl0);

  public:
  virtual void
      evaluate(seissol::initializers::Layer& layerData,
               seissol::initializers::DynamicRupture* dynRup,
               real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
               real fullUpdateTime,
               double timeWeights[CONVERGENCE_ORDER]) override;
};

#endif // SEISSOL_AGINGLAW_H

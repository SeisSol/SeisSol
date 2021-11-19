#ifndef SEISSOL_AGINGLAW_H
#define SEISSOL_AGINGLAW_H

#include "BaseFrictionLaw.h"

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class AgingLaw : public BaseFrictionLaw<AgingLaw> {
  public:
  using BaseFrictionLaw::BaseFrictionLaw;
  using BaseFrictionLaw::copyLtsTreeToLocal;

  protected:
  real calcStateVariableHook(real SV0, real tmp, real time_inc, real rs_sl0);

  public:
  void evaluate(seissol::initializers::Layer& layerData,
                seissol::initializers::DynamicRupture* dynRup,
                real (*QInterpolatedPlus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real (*QInterpolatedMinus)[CONVERGENCE_ORDER][tensor::QInterpolated::size()],
                real fullUpdateTime,
                double timeWeights[CONVERGENCE_ORDER]) override;
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_AGINGLAW_H

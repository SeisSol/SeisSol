#ifndef SEISSOL_SLIPLAW_H
#define SEISSOL_SLIPLAW_H

#include "SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL4. Since FL4 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw>::copyLtsTreeToLocal;

  /**
   * Integrates the state variable ODE in time
   * \f[\frac{\partial \Theta}{\partial t} = - \frac{V}{L}\Theta \cdot \log\left( \frac{V}{L} \Theta
   * \right). \f]
   * Analytic solution
   * \f[ \Theta(t) = \frac{L}{V} \exp\left[ \log\left( \frac{V
   * \Theta_0}{L}\right) \exp\left( - \frac{V}{L} t\right)\right].\f]
   * @param stateVarReference \f$ \Theta_0 \f$
   * @param timeIncremetn \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Theta(t) \f$
   */
  real updateStateVariable(int pointIndex,
                           unsigned int face,
                           real stateVarReference,
                           real timeIncrement,
                           real localSlipRate);
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_SLIPLAW_H

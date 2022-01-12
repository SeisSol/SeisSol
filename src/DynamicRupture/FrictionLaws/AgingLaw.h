#ifndef SEISSOL_AGINGLAW_H
#define SEISSOL_AGINGLAW_H

#include "SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
class AgingLaw : public SlowVelocityWeakeningLaw<AgingLaw> {
  public:
  using SlowVelocityWeakeningLaw<AgingLaw>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<AgingLaw>::copyLtsTreeToLocal;

  /**
   * Integrates the state variable ODE in time
   * \f[ \frac{\partial \Theta}{\partial t} = 1 - \frac{V}{L} \Theta. \f]
   * Analytic solution:
   * \f[\Theta(t) = - \Theta_0 \frac{V}{L} \cdot \exp\left( -\frac{V}{L} \cdot t\right) + \exp\left(
   * -\frac{V}{L} \cdot t\right). \f]
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
#endif // SEISSOL_AGINGLAW_H

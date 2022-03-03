#ifndef SEISSOL_SLIPLAW_H
#define SEISSOL_SLIPLAW_H

#include "SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law {

/**
 * This class was not tested and compared to the Fortran FL4. Since FL4 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
template <typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::copyLtsTreeToLocal;

  /**
   * Integrates the state variable ODE in time
   * \f[\frac{\partial \Theta}{\partial t} = - \frac{V}{L}\Theta \cdot \log\left( \frac{V}{L} \Theta
   * \right). \f]
   * Analytic solution
   * \f[ \Theta(t) = \frac{L}{V} \exp\left[ \log\left( \frac{V
   * \Theta_0}{L}\right) \exp\left( - \frac{V}{L} t\right)\right].\f]
   * Note that we need double precision here, since single precision led to NaNs.
   * @param stateVarReference \f$ \Theta_0 \f$
   * @param timeIncremetn \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Theta(t) \f$
   */
  double updateStateVariable(int pointIndex,
                             unsigned int face,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    double localSl0 = this->sl0[face][pointIndex];
    double exp1 = exp(-localSlipRate * (timeIncrement / localSl0));
    return localSl0 / localSlipRate * std::pow(localSlipRate * stateVarReference / localSl0, exp1);
  }
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_SLIPLAW_H

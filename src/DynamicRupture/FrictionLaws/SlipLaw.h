#ifndef SEISSOL_SLIPLAW_H
#define SEISSOL_SLIPLAW_H

#include "SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law {
template <typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::copyLtsTreeToLocal;

  /**
   * Integrates the state variable ODE in time
   * \f[\frac{\partial \Psi}{\partial t} = - \frac{V}{L}\Psi \cdot \log\left( \frac{V}{L} \Psi
   * \right). \f]
   * Analytic solution
   * \f[ \Psi(t) = \frac{L}{V} \exp\left[ \log\left( \frac{V
   * \Psi_0}{L}\right) \exp\left( - \frac{V}{L} t\right)\right].\f]
   * Note that we need double precision here, since single precision led to NaNs.
   * @param stateVarReference \f$ \Psi_0 \f$
   * @param timeIncremetn \f$ t \f$
   * @param localSlipRate \f$ V \f$
   * @return \f$ \Psi(t) \f$
   */
  #pragma omp declare simd
  double updateStateVariable(int pointIndex,
                             unsigned int face,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    const double localSl0 = this->sl0[face][pointIndex];
    const double exp1 = std::exp(-localSlipRate * (timeIncrement / localSl0));
    return localSl0 / localSlipRate * std::pow(localSlipRate * stateVarReference / localSl0, exp1);
  }
};

} // namespace seissol::dr::friction_law
#endif // SEISSOL_SLIPLAW_H

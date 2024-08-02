#ifndef SEISSOL_ADJOINTSLIP_H
#define SEISSOL_ADJOINTSLIP_H

#include "AdjointSlowVelWeakening.h"

namespace seissol::dr::friction_law {

/**
 * This is inherited from SlipLaws, including the following descriptions.
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
template <class TPMethod>
class AdjointSlip : public AdjointSlowVelWeakening<AdjointSlip<TPMethod>, TPMethod> {
  public:
  using AdjointSlowVelWeakening<AdjointSlip<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using AdjointSlowVelWeakening<AdjointSlip<TPMethod>, TPMethod>::copyLtsTreeToLocal;

/**
 * Integrates the state variable ODE in time
 * \f[ \frac{\partial \Psi}{\partial t} = 1 - \frac{V}{L} \Psi \f]
 * Analytic solution:
 * \f[\Psi(t) = - \Psi_0 \frac{V}{L} \cdot \exp\left( -\frac{V}{L} \cdot t\right) + \exp\left(
 * -\frac{V}{L} \cdot t\right). \f]
 * Note that we need double precision here, since single precision led to NaNs.
 * @param stateVarReference \f$ \Psi_0 \f$
 * @param timeIncrement \f$ t \f$
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
#endif // SEISSOL_ADJOINTSLIP_H

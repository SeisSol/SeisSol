// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_AGINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_AGINGLAW_H_

#include "SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::cpu {

/**
 * This class was not tested and compared to the Fortran FL3. Since FL3 initialization did not work
 * properly on the Master Branch. This class is also less optimized. It was left in here to have a
 * reference of how it could be implemented.
 */
template <class TPMethod>
class AgingLaw : public SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod>::copyLtsTreeToLocal;

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
  [[nodiscard]] double updateStateVariable(int pointIndex,
                                           unsigned int face,
                                           double stateVarReference,
                                           double timeIncrement,
                                           double localSlipRate) const {
    const double localSl0 = this->sl0[face][pointIndex];
    const double preexp1 = -localSlipRate * (timeIncrement / localSl0);
    const double exp1 = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    return stateVarReference * exp1 + localSl0 / localSlipRate * exp1m;
  }
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_AGINGLAW_H_

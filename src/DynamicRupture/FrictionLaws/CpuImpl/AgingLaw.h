// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_AGINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_AGINGLAW_H_

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
  using SlowVelocityWeakeningLaw<AgingLaw<TPMethod>, TPMethod>::copyStorageToLocal;

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
  /// generic over the scalar the slip rate arrives in; see SlowVelocityWeakeningLaw::StateScalar
  template <typename S>
  [[nodiscard]] S updateStateVariable(std::uint32_t pointIndex,
                                      std::size_t faceIndex,
                                      double stateVarReference,
                                      double timeIncrement,
                                      S localSlipRate) const {
    using std::exp;
    using std::expm1;
    const double localSl0 = this->sl0_[faceIndex][pointIndex];
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    const S exp1v = exp(preexp1);
    const S exp1m = -expm1(preexp1);
    return S(stateVarReference) * exp1v + S(localSl0) / localSlipRate * exp1m;
  }
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_AGINGLAW_H_

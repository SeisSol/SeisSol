// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLIPLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLIPLAW_H_

#include "SlowVelocityWeakeningLaw.h"

namespace seissol::dr::friction_law::cpu {
template <typename TPMethod>
class SlipLaw : public SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod> {
  public:
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::SlowVelocityWeakeningLaw;
  using SlowVelocityWeakeningLaw<SlipLaw<TPMethod>, TPMethod>::copyStorageToLocal;

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
  /// generic over the scalar the slip rate arrives in; see SlowVelocityWeakeningLaw::StateScalar
  template <typename S>
  S updateStateVariable(std::uint32_t pointIndex,
                        std::size_t faceIndex,
                        real stateVarReference,
                        real timeIncrement,
                        S localSlipRate) {
    using std::exp;
    using std::expm1;
    using std::log;
    const real localSl0 = this->sl0_[faceIndex][pointIndex];
    const S preexp1 = -localSlipRate * S(timeIncrement / localSl0);
    const S exp1v = exp(preexp1);
    const S exp1m = -expm1(preexp1);
    // (L / V) (V Psi / L)^e is the weighted geometric mean of L / V and Psi, with the weights
    // 1 - e and e, so it is an exponential of logarithms and nothing else. Formed that way the
    // large quotient and the small power beside it never meet. As a product they have to: each
    // factor's derivative is of the order L / V^2 while their product's is of the order t / V,
    // fourteen decades below it at a slip rate of a picometre a second, so in single precision the
    // derivative comes back as exactly zero for every slip rate under a millimetre a second and as
    // an infinity under 1e-22. Both the base and the exponent follow the slip rate, and here both
    // do so through the sum in the exponent.
    //
    // The first weight goes through expm1, because 1 - exp(-z) cancels wherever the relaxation is
    // slight and the logarithm it multiplies reaches eighty: the error would land undivided in the
    // exponent, five parts in a thousand of the state at a slip rate of a micrometre a second.
    const S logQuotient = S(std::log(localSl0)) - log(localSlipRate);
    return exp(exp1m * logQuotient + exp1v * S(std::log(stateVarReference)));
  }
};

} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLIPLAW_H_

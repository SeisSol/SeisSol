// SPDX-FileCopyrightText: 2022-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SEVEREVELOCITYWEAKENINGLAW_H_

#include "RateAndState.h"

namespace seissol::dr::friction_law::cpu {
template <class TPMethod>
class SevereVelocityWeakeningLaw
    : public RateAndStateBase<SevereVelocityWeakeningLaw<TPMethod>, TPMethod> {
  public:
  using RateAndStateBase<SevereVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyLtsTreeToLocal(seissol::initializer::Layer& layerData,
                          const seissol::initializer::DynamicRupture* const dynRup,
                          real fullUpdateTime) {}

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  double updateStateVariable(int pointIndex,
                             unsigned int face,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    const double localSl0 = this->sl0[face][pointIndex];

    const double steadyStateStateVariable = localSlipRate * localSl0 / this->drParameters->rsSr0;

    const double preexp1 = -this->drParameters->rsSr0 * (timeIncrement / localSl0);
    const double exp1 = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const double localStateVariable = steadyStateStateVariable * exp1m + exp1 * stateVarReference;

    return localStateVariable;
  }

/**
 * Computes the friction coefficient from the state variable and slip rate
 * \f[\mu = a \cdot \sinh^{-1} \left( \frac{V}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi
 * / L)}{a} \right)\right).\f]
 * Note that we need double precision here, since single precision led to NaNs.
 * @param localSlipRateMagnitude \f$ V \f$
 * @param localStateVariable \f$ \Psi \f$
 * @return \f$ \mu \f$
 */
#pragma omp declare simd
  double updateMu(unsigned int ltsFace,
                  unsigned int pointIndex,
                  double localSlipRateMagnitude,
                  double localStateVariable) {
    const double localA = this->a[ltsFace][pointIndex];
    const double localSl0 = this->sl0[ltsFace][pointIndex];
    const double c = this->drParameters->rsB * localStateVariable / (localStateVariable + localSl0);
    return this->drParameters->rsF0 +
           localA * localSlipRateMagnitude / (localSlipRateMagnitude + this->drParameters->rsSr0) -
           c;
  }

/**
 * Computes the derivative of the friction coefficient with respect to the slip rate.
 * \f[\frac{\partial}{\partial V}\mu = \frac{aC}{\sqrt{(VC)^2 +1}} \text{ with } C =
 * \frac{1}{2V_0} \cdot \exp \left(\frac{f_0 + b \log(V_0\Psi / L)}{a} \right). \f]
 * Note that we need double precision here, since single precision led to NaNs.
 * @param localSlipRateMagnitude \f$ V \f$
 * @param localStateVariable \f$ \Psi \f$
 * @return \f$ \mu \f$
 */
#pragma omp declare simd
  double updateMuDerivative(unsigned int ltsFace,
                            unsigned int pointIndex,
                            double localSlipRateMagnitude,
                            double localStateVariable) {
    const double localA = this->a[ltsFace][pointIndex];
    const double divisor = (localSlipRateMagnitude + this->drParameters->rsSr0);
    return localA * this->drParameters->rsSr0 / (divisor * divisor);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                        unsigned int ltsFace) const {
#pragma omp simd
    for (size_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] = stateVariableBuffer[pointIndex];
    }
  }

  void executeIfNotConverged(const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                             unsigned ltsFace) {}
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_SEVEREVELOCITYWEAKENINGLAW_H_

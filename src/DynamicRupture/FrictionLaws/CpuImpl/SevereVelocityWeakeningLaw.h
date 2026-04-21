// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

#include "DynamicRupture/Misc.h"
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
  void copyStorageToLocal(DynamicRupture::Layer& layerData) {}

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  double updateStateVariable(std::uint32_t pointIndex,
                             std::size_t faceIndex,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    const real localSl0 = this->sl0_[faceIndex][pointIndex];

    const real steadyStateStateVariable = localSlipRate * localSl0 / this->drParameters_->rsSr0;

    const double preexp1 = -this->drParameters_->rsSr0 * (timeIncrement / localSl0);
    const double exp1v = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const double localStateVariable = steadyStateStateVariable * exp1m + exp1v * stateVarReference;

    return localStateVariable;
  }

  struct MuDetails {
    std::array<real, misc::NumPaddedPoints> a{};
    std::array<real, misc::NumPaddedPoints> c{};
    std::array<real, misc::NumPaddedPoints> f0{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const real localA = this->a_[ltsFace][pointIndex];
      const real localSl0 = this->sl0_[ltsFace][pointIndex];
      const real c = this->b_[ltsFace][pointIndex] * localStateVariable[pointIndex] /
                     (localStateVariable[pointIndex] + localSl0);

      details.a[pointIndex] = localA;
      details.c[pointIndex] = c;
      details.f0[pointIndex] = this->f0_[ltsFace][pointIndex];
    }
    return details;
  }

#pragma omp declare simd
  real updateMu(std::uint32_t pointIndex, real localSlipRateMagnitude, const MuDetails& details) {
    return details.f0[pointIndex] +
           details.a[pointIndex] * localSlipRateMagnitude /
               (localSlipRateMagnitude + this->drParameters_->rsSr0) -
           details.c[pointIndex];
  }

#pragma omp declare simd
  real updateMuDerivative(std::uint32_t pointIndex,
                          real localSlipRateMagnitude,
                          const MuDetails& details) {
    // note that: d/dx (x/(x+c)) = ((x+c)-x)/(x+c)**2 = c/(x+c)**2
    const real divisor = (localSlipRateMagnitude + this->drParameters_->rsSr0);
    return details.a[pointIndex] * this->drParameters_->rsSr0 / (divisor * divisor);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                        std::size_t ltsFace) const {
#pragma omp simd
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->stateVariable_[ltsFace][pointIndex] = stateVariableBuffer[pointIndex];
    }
  }

  void executeIfNotConverged(const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                             std::size_t ltsFace) {}
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

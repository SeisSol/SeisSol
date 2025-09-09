// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

#include "RateAndState.h"
#include <DynamicRupture/Misc.h>

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
                          const seissol::initializer::DynamicRupture* const dynRup) {}

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  double updateStateVariable(std::uint32_t pointIndex,
                             std::size_t faceIndex,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    const double localSl0 = this->sl0[faceIndex][pointIndex];

    const double steadyStateStateVariable = localSlipRate * localSl0 / this->drParameters->rsSr0;

    const double preexp1 = -this->drParameters->rsSr0 * (timeIncrement / localSl0);
    const double exp1v = std::exp(preexp1);
    const double exp1m = -std::expm1(preexp1);
    const double localStateVariable = steadyStateStateVariable * exp1m + exp1v * stateVarReference;

    return localStateVariable;
  }

  struct MuDetails {
    std::array<double, misc::NumPaddedPoints> a{};
    std::array<double, misc::NumPaddedPoints> c{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; ++pointIndex) {
      const double localA = this->a[ltsFace][pointIndex];
      const double localSl0 = this->sl0[ltsFace][pointIndex];
      const double c = this->drParameters->rsB *
                       static_cast<double>(localStateVariable[pointIndex]) /
                       (static_cast<double>(localStateVariable[pointIndex]) + localSl0);
      details.a[pointIndex] = localA;
      details.c[pointIndex] = c;
    }
    return details;
  }

#pragma omp declare simd
  double
      updateMu(std::uint32_t pointIndex, double localSlipRateMagnitude, const MuDetails& details) {
    return this->drParameters->rsF0 +
           details.a[pointIndex] * localSlipRateMagnitude /
               (localSlipRateMagnitude + this->drParameters->rsSr0) -
           details.c[pointIndex];
  }

#pragma omp declare simd
  double updateMuDerivative(std::uint32_t pointIndex,
                            double localSlipRateMagnitude,
                            const MuDetails& details) {
    // note that: d/dx (x/(x+c)) = ((x+c)-x)/(x+c)**2 = c/(x+c)**2
    const double divisor = (localSlipRateMagnitude + this->drParameters->rsSr0);
    return details.a[pointIndex] * this->drParameters->rsSr0 / (divisor * divisor);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints>& stateVariableBuffer,
                        std::size_t ltsFace) const {
#pragma omp simd
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] = stateVariableBuffer[pointIndex];
    }
  }

  void executeIfNotConverged(const std::array<real, misc::NumPaddedPoints>& localStateVariable,
                             std::size_t ltsFace) {}
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SEVEREVELOCITYWEAKENINGLAW_H_

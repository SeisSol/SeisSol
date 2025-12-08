// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_

#include "RateAndState.h"

namespace seissol::dr::friction_law::cpu {
template <typename Cfg, class Derived, class TPMethod>
class SlowVelocityWeakeningLaw
    : public RateAndStateBase<Cfg, SlowVelocityWeakeningLaw<Cfg, Derived, TPMethod>, TPMethod> {
  public:
  using real = Real<Cfg>;
  SlowVelocityWeakeningLaw() = default;
  using RateAndStateBase<Cfg, SlowVelocityWeakeningLaw, TPMethod>::RateAndStateBase;

  /**
   * copies all parameters from the DynamicRupture LTS to the local attributes
   */
  void copyStorageToLocal(DynamicRupture::Layer& layerData) override {}

  std::unique_ptr<FrictionSolver> clone() override {
    return std::make_unique<Derived>(*static_cast<Derived*>(this));
  }

// Note that we need double precision here, since single precision led to NaNs.
#pragma omp declare simd
  double updateStateVariable(std::uint32_t pointIndex,
                             std::size_t faceIndex,
                             double stateVarReference,
                             double timeIncrement,
                             double localSlipRate) {
    return static_cast<Derived*>(this)->updateStateVariable(
        pointIndex, faceIndex, stateVarReference, timeIncrement, localSlipRate);
  }

  struct MuDetails {
    std::array<real, misc::NumPaddedPoints<Cfg>> a{};
    std::array<real, misc::NumPaddedPoints<Cfg>> cLin{};
    std::array<real, misc::NumPaddedPoints<Cfg>> cExpLog{};
    std::array<real, misc::NumPaddedPoints<Cfg>> cExp{};
    std::array<real, misc::NumPaddedPoints<Cfg>> acLin{};
  };

  MuDetails getMuDetails(std::size_t ltsFace,
                         const std::array<real, misc::NumPaddedPoints<Cfg>>& localStateVariable) {
    MuDetails details{};
#pragma omp simd
    for (std::uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; ++pointIndex) {
      const real localA = this->a[ltsFace][pointIndex];
      const real localSl0 = this->sl0[ltsFace][pointIndex];
      const real log1 =
          std::log(this->drParameters->rsSr0 * localStateVariable[pointIndex] / localSl0);
      const real localF0 = this->f0[ltsFace][pointIndex];
      const real localB = this->b[ltsFace][pointIndex];

      const real cLin = 0.5 / this->drParameters->rsSr0;
      const real cExpLog = (localF0 + localB * log1) / localA;
      const real cExp = rs::computeCExp(cExpLog);
      const real acLin = localA * cLin;

      details.a[pointIndex] = localA;
      details.cLin[pointIndex] = cLin;
      details.cExpLog[pointIndex] = cExpLog;
      details.cExp[pointIndex] = cExp;
      details.acLin[pointIndex] = acLin;
    }
    return details;
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
  real updateMu(std::uint32_t pointIndex, real localSlipRateMagnitude, const MuDetails& details) {
    const real lx = details.cLin[pointIndex] * localSlipRateMagnitude;
    return details.a[pointIndex] *
           rs::arsinhexp(lx, details.cExpLog[pointIndex], details.cExp[pointIndex]);
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
  real updateMuDerivative(std::uint32_t pointIndex,
                          real localSlipRateMagnitude,
                          const MuDetails& details) {
    const real lx = details.cLin[pointIndex] * localSlipRateMagnitude;
    return details.acLin[pointIndex] *
           rs::arsinhexpDerivative(lx, details.cExpLog[pointIndex], details.cExp[pointIndex]);
  }

  /**
   * Resample the state variable. For Slow Velocity Weakening Laws, we just copy the buffer into the
   * member variable.
   */
  void resampleStateVar(const std::array<real, misc::NumPaddedPoints<Cfg>>& stateVariableBuffer,
                        std::size_t ltsFace) const {
#pragma omp simd
    for (uint32_t pointIndex = 0; pointIndex < misc::NumPaddedPoints<Cfg>; pointIndex++) {
      this->stateVariable[ltsFace][pointIndex] = stateVariableBuffer[pointIndex];
    }
  }

  void executeIfNotConverged(const std::array<real, misc::NumPaddedPoints<Cfg>>& localStateVariable,
                             std::size_t ltsFace) {
    [[maybe_unused]] const real tmp =
        0.5 / this->drParameters->rsSr0 *
        std::exp(
            (this->drParameters->rsF0 +
             this->drParameters->rsB * std::log(this->drParameters->rsSr0 * localStateVariable[0] /
                                                this->drParameters->rsSr0)) /
            this->a[ltsFace][0]);
    assert(!std::isnan(tmp) && "nonConvergence RS Newton");
  }
};
} // namespace seissol::dr::friction_law::cpu

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FRICTIONLAWS_CPUIMPL_SLOWVELOCITYWEAKENINGLAW_H_

#include "LinearSlipWeakening.h"
namespace seissol::dr::friction_law {

void ForcedRuptureTime::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                           seissol::initializers::DynamicRupture* dynRup,
                                           real fullUpdateTime) {
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime*>(dynRup);
  forcedRuptureTime = layerData.var(concreteLts->forcedRuptureTime);
}

void ForcedRuptureTime::stateVariableHook(std::array<real, misc::numPaddedPoints>& stateVariable,
                                          real time,
                                          unsigned int ltsFace) {
  for (unsigned pointIndex = 0; pointIndex < misc::numPaddedPoints; pointIndex++) {
    real f2 = 0.0;
    if (drParameters.t0 == 0) {
      if (time >= forcedRuptureTime[ltsFace][pointIndex]) {
        f2 = 1.0;
      } else {
        f2 = 0.0;
      }
    } else {
      f2 = std::max(static_cast<real>(0.0),
                    std::min(static_cast<real>(1.0),
                             // Note: In the fortran implementation on the master branch, this is
                             // m_fullUpdateTime, but this implementation is correct.
                             (time - forcedRuptureTime[ltsFace][pointIndex]) / drParameters.t0));
    }
    stateVariable[pointIndex] = std::max(stateVariable[pointIndex], f2);
  }
}

void BiMaterialFault::copyLtsTreeToLocal(seissol::initializers::Layer& layerData,
                                         seissol::initializers::DynamicRupture* dynRup,
                                         real fullUpdateTime) {
  // maybe change later to const_cast?
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTS_LinearSlipWeakeningBimaterial*>(dynRup);
  regularisedStrength = layerData.var(concreteLts->regularisedStrength);
}

void BiMaterialFault::strengthHook(real& strength,
                                   real& localSlipRate,
                                   real& sigma,
                                   real& mu,
                                   real& deltaT,
                                   unsigned int ltsFace,
                                   unsigned int pointIndex) {
  // modify strength according to prakash clifton
  // see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  prak_clif_mod(strength, sigma, localSlipRate, mu, deltaT);
  // save for output
  regularisedStrength[ltsFace][pointIndex] = strength;
}

/*
 * calculates prakash clifton regularization
 */
void BiMaterialFault::prak_clif_mod(
    real& strength, real& sigma, real& localSlipRate, real& mu, real& dt) {
  real expterm =
      std::exp(-(std::abs(localSlipRate) + drParameters.vStar) * dt / drParameters.prakashLength);
  strength = strength * expterm - std::max(static_cast<real>(0.0), -mu * sigma) * (expterm - 1.0);
}
} // namespace seissol::dr::friction_law
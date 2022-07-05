#include "LinearSlipWeakening.h"
namespace seissol::dr::friction_law {

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
  // modify strength according to Prakash-Clifton
  // see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  prak_clif_mod(strength, sigma, localSlipRate, mu, deltaT);
  // save for output
  regularisedStrength[ltsFace][pointIndex] = strength;
}

/*
 * calculates Prakash-Clifton regularization
 */
void BiMaterialFault::prak_clif_mod(
    real& strength, real& sigma, real& localSlipRate, real& mu, real& dt) {
  real expterm =
      std::exp(-(std::abs(localSlipRate) + drParameters.vStar) * dt / drParameters.prakashLength);
  strength = strength * expterm - std::max(static_cast<real>(0.0), -mu * sigma) * (expterm - 1.0);
}
} // namespace seissol::dr::friction_law
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

real BiMaterialFault::strengthHook(real& faultStrength,
                                   real& localSlipRate,
                                   real& deltaT,
                                   unsigned int ltsFace,
                                   unsigned int pointIndex) {
  // modify strength according to Prakash-Clifton
  // see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  real newStrength =
      prak_clif_mod(regularisedStrength[ltsFace][pointIndex], faultStrength, localSlipRate, deltaT);
  regularisedStrength[ltsFace][pointIndex] = newStrength;

  return newStrength;
}

real BiMaterialFault::prak_clif_mod(real& regularisedStrength,
                                    real& faultStrength,
                                    real& localSlipRate,
                                    real& dt) {
  real expterm =
      std::exp(-(std::abs(localSlipRate) + drParameters.vStar) * dt / drParameters.prakashLength);
  real newstrength = regularisedStrength * expterm + faultStrength * (1.0 - expterm);
  return newstrength;
}
} // namespace seissol::dr::friction_law
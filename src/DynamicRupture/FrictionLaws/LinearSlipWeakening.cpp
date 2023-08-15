#include "LinearSlipWeakening.h"
#include "Common/configtensor.hpp"
namespace seissol::dr::friction_law {

template <typename Config>
void NoSpecialization<Config>::resampleSlipRate(
    RealT (&resampledSlipRate)[dr::misc::numPaddedPoints<Config>],
    RealT const (&slipRateMagnitude)[dr::misc::numPaddedPoints<Config>]) {
  typename Yateto<Config>::Kernel::dynamicRupture::resampleParameter resampleKrnl;
  resampleKrnl.resample = Yateto<Config>::Init::resample::Values;
  resampleKrnl.originalQ = slipRateMagnitude;
  resampleKrnl.resampledQ = resampledSlipRate;
  resampleKrnl.execute();
}
template <typename Config>
void BiMaterialFault<Config>::copyLtsTreeToLocal(
    seissol::initializers::Layer& layerData,
    seissol::initializers::DynamicRupture<Config> const* const dynRup,
    RealT fullUpdateTime) {
  auto* concreteLts =
      dynamic_cast<seissol::initializers::LTSLinearSlipWeakeningBimaterial<Config> const* const>(
          dynRup);
  regularisedStrength = layerData.var(concreteLts->regularisedStrength);
}

#pragma omp declare simd
template <typename Config>
typename Config::RealT BiMaterialFault<Config>::strengthHook(RealT faultStrength,
                                                             RealT localSlipRate,
                                                             RealT deltaT,
                                                             unsigned int ltsFace,
                                                             unsigned int pointIndex) {
  // modify strength according to Prakash-Clifton
  // see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  const RealT expterm =
      std::exp(-(std::max(static_cast<RealT>(0.0), localSlipRate) + drParameters->vStar) * deltaT /
               drParameters->prakashLength);
  const RealT newStrength =
      regularisedStrength[ltsFace][pointIndex] * expterm + faultStrength * (1.0 - expterm);
  regularisedStrength[ltsFace][pointIndex] = newStrength;
  return newStrength;
}

} // namespace seissol::dr::friction_law

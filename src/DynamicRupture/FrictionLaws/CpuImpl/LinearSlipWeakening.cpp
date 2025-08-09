// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LinearSlipWeakening.h"
#include "DynamicRupture/Misc.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Layer.h"
#include <GeneratedCode/init.h>
#include <GeneratedCode/kernel.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace seissol::dr::friction_law::cpu {

template <typename Cfg>
void NoSpecialization<Cfg>::resampleSlipRate(
    real (&resampledSlipRate)[dr::misc::NumPaddedPoints<Cfg>],
    const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints<Cfg>]) {
  dynamicRupture::kernel::resampleParameter<Cfg> resampleKrnl;
  resampleKrnl.resample = init::resample<Cfg>::Values;
  resampleKrnl.originalQ = slipRateMagnitude;
  resampleKrnl.resampledQ = resampledSlipRate;
  resampleKrnl.execute();
}

template <typename Cfg>
void BiMaterialFault<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  regularizedStrength = layerData.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(Cfg());
}

#pragma omp declare simd
template <typename Cfg>
Real<Cfg> BiMaterialFault<Cfg>::strengthHook(real faultStrength,
                                             real localSlipRate,
                                             real deltaT,
                                             std::size_t ltsFace,
                                             std::uint32_t pointIndex) {
  // modify strength according to Prakash-Clifton
  // see e.g.: Pelties - Verification of an ADER-DG method for complex dynamic rupture problems
  const real expterm =
      std::exp(-(std::max(static_cast<real>(0.0), localSlipRate) + drParameters->vStar) * deltaT /
               drParameters->prakashLength);
  const real newStrength =
      regularizedStrength[ltsFace][pointIndex] * expterm + faultStrength * (1.0 - expterm);
  regularizedStrength[ltsFace][pointIndex] = newStrength;
  return newStrength;
}

#pragma omp declare simd
template <typename Cfg>
Real<Cfg> TPApprox<Cfg>::stateVariableHook(real localAccumulatedSlip,
                                           real localDc,
                                           std::size_t ltsFace,
                                           std::uint32_t pointIndex) {
  const real factor = (1.0 + std::fabs(localAccumulatedSlip) / localDc);
  return 1.0 - std::pow(factor, -drParameters->tpProxyExponent);
}

#define SEISSOL_CONFIGITER(cfg) template class NoSpecialization<cfg>;
#include "ConfigInclude.h"
#define SEISSOL_CONFIGITER(cfg) template class BiMaterialFault<cfg>;
#include "ConfigInclude.h"
#define SEISSOL_CONFIGITER(cfg) template class TPApprox<cfg>;
#include "ConfigInclude.h"

} // namespace seissol::dr::friction_law::cpu

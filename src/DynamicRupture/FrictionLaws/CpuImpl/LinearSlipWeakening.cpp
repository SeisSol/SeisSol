// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LinearSlipWeakening.h"

#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::cpu {

template <typename Cfg>
void NoSpecialization<Cfg>::allocateAuxiliaryMemory(GlobalData<Cfg>* globalData) {
  resampleKrnlPrototype_.bindGlobals(*globalData);
}

template <typename Cfg>
void NoSpecialization<Cfg>::resampleSlipRate(
    real (&resampledSlipRate)[dr::misc::NumPaddedPoints<Cfg>],
    const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints<Cfg>]) const {
  auto resampleKrnl = resampleKrnlPrototype_;
  resampleKrnl.originalQ = slipRateMagnitude;
  resampleKrnl.resampledQ = resampledSlipRate;
  resampleKrnl.execute();
}
template <typename Cfg>
void BiMaterialFault<Cfg>::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  regularizedStrength_ =
      layerData.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(Cfg());
}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template class NoSpecialization<Cfg>;                                                            \
  template class BiMaterialFault<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::dr::friction_law::cpu

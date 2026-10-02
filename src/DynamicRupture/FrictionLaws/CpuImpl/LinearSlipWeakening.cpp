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
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::cpu {

void NoSpecialization::allocateAuxiliaryMemory(GlobalData<Config>* globalData) {
  resampleKrnlPrototype_.bindGlobals(*globalData);
}

void NoSpecialization::resampleSlipRate(
    real (&resampledSlipRate)[dr::misc::NumPaddedPoints<Config>],
    const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints<Config>]) const {
  auto resampleKrnl = resampleKrnlPrototype_;
  resampleKrnl.originalQ = slipRateMagnitude;
  resampleKrnl.resampledQ = resampledSlipRate;
  resampleKrnl.execute();
}
void BiMaterialFault::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  regularizedStrength_ =
      layerData.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>(Config());
}

} // namespace seissol::dr::friction_law::cpu

// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "LinearSlipWeakening.h"

#include "DynamicRupture/Misc.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"

namespace seissol::dr::friction_law::cpu {

void NoSpecialization::resampleSlipRate(
    real (&resampledSlipRate)[dr::misc::NumPaddedPoints],
    const real (&slipRateMagnitude)[dr::misc::NumPaddedPoints]) {
  dynamicRupture::kernel::resampleParameter resampleKrnl;
  resampleKrnl.resample = init::resample::Values;
  resampleKrnl.originalQ = slipRateMagnitude;
  resampleKrnl.resampledQ = resampledSlipRate;
  resampleKrnl.execute();
}
void BiMaterialFault::copyStorageToLocal(DynamicRupture::Layer& layerData) {
  regularizedStrength_ = layerData.var<LTSLinearSlipWeakeningBimaterial::RegularizedStrength>();
}

} // namespace seissol::dr::friction_law::cpu

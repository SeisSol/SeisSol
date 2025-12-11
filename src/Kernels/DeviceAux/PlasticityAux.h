// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

#include "Equations/Datastructures.h"
#include "Initializer/BasicTypedefs.h"
#include "Model/Plasticity.h"

#include <stddef.h>

namespace seissol::kernels::device::aux::plasticity {
constexpr static int NumStressComponents = model::MaterialT::TractionQuantities;

void plasticityNonlinear(real** __restrict nodalStressTensors,
                         real** __restrict prevNodal,
                         real** __restrict pstrainPtr,
                         unsigned* __restrict isAdjustableVector,
                         std::size_t* __restrict yieldCounter,
                         const seissol::model::PlasticityData* __restrict plasticity,
                         double oneMinusIntegratingFactor,
                         double tV,
                         double timeStepWidth,
                         size_t numElements,
                         void* streamPtr);
} // namespace seissol::kernels::device::aux::plasticity

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

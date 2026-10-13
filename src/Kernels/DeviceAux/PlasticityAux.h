// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

#include "Config.h"
#include "Equations/Datastructures.h"
#include "Initializer/BasicTypedefs.h"
#include "Model/Plasticity.h"

#include <stddef.h>

namespace seissol::kernels::device::aux::plasticity {
/// The stress components of the material of the configuration `Cfg`, which plasticity adjusts.
template <typename Cfg>
constexpr int NumStressComponents = model::MaterialOf<Cfg>::TractionComponents;

template <typename Cfg>
void plasticityNonlinear(Real<Cfg>** __restrict nodalStressTensors,
                         Real<Cfg>** __restrict pstrainPtr,
                         unsigned* __restrict isAdjustableVector,
                         std::size_t* __restrict yieldCounter,
                         const seissol::model::PlasticityData<Cfg>* __restrict plasticity,
                         Real<Cfg> oneMinusIntegratingFactor,
                         Real<Cfg> tV,
                         Real<Cfg> timeStepWidth,
                         size_t numElements,
                         void* streamPtr);
} // namespace seissol::kernels::device::aux::plasticity

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

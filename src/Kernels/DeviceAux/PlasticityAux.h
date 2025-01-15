// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_
#define SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

#include "Initializer/BasicTypedefs.h"
#include "Model/Plasticity.h"
#include <stddef.h>

#define NUM_STRESS_COMPONENTS 6

// NOTE: using c++14 because of cuda@10
namespace seissol {
namespace kernels {
namespace device {
namespace aux {
namespace plasticity {
void adjustDeviatoricTensors(real** nodalStressTensors,
                             unsigned* isAdjustableVector,
                             const seissol::model::PlasticityData* plasticity,
                             double oneMinusIntegratingFactor,
                             size_t numElements,
                             void* streamPtr);

void adjustPointers(real* qEtaNodal,
                    real** qEtaNodalPtrs,
                    real* qEtaModal,
                    real** qEtaModalPtrs,
                    real* dUdTpstrain,
                    real** dUdTpstrainPtrs,
                    size_t numElements,
                    void* streamPtr);

void computePstrains(real** pstrains,
                     const seissol::model::PlasticityData* plasticityData,
                     real** dofs,
                     real* prevDofs,
                     real** dUdTpstrain,
                     double tV,
                     double oneMinusIntegratingFactor,
                     double timeStepWidth,
                     unsigned* isAdjustableVector,
                     size_t numElements,
                     void* streamPtr);

void pstrainToQEtaModal(real** pstrains,
                        real** qEtaModalPtrs,
                        unsigned* isAdjustableVector,
                        size_t numElements,
                        void* streamPtr);

void qEtaModalToPstrain(real** qEtaModalPtrs,
                        real** pstrains,
                        unsigned* isAdjustableVector,
                        size_t numElements,
                        void* streamPtr);

void updateQEtaNodal(real** qEtaNodalPtrs,
                     real** qStressNodalPtrs,
                     double timeStepWidth,
                     unsigned* isAdjustableVector,
                     size_t numElements,
                     void* streamPtr);

} // namespace plasticity
} // namespace aux
} // namespace device
} // namespace kernels
} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

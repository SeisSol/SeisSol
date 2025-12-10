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

namespace seissol::kernels::device::aux {

template <typename Cfg>
struct Plasticity {
  using real = Real<Cfg>;

  static void adjustDeviatoricTensors(real** nodalStressTensors,
                                      unsigned* isAdjustableVector,
                                      const seissol::model::PlasticityData<Cfg>* plasticity,
                                      double oneMinusIntegratingFactor,
                                      size_t numElements,
                                      void* streamPtr);

  static void computePstrains(real** pstrains,
                              const seissol::model::PlasticityData<Cfg>* plasticityData,
                              real** dofs,
                              real** prevDofs,
                              real** dUdTpstrain,
                              double tV,
                              double oneMinusIntegratingFactor,
                              double timeStepWidth,
                              unsigned* isAdjustableVector,
                              size_t numElements,
                              void* streamPtr);

  static void updateQEtaNodal(real** qEtaNodalPtrs,
                              real** qStressNodalPtrs,
                              double timeStepWidth,
                              unsigned* isAdjustableVector,
                              size_t numElements,
                              void* streamPtr);
};
} // namespace seissol::kernels::device::aux

#endif // SEISSOL_SRC_KERNELS_DEVICEAUX_PLASTICITYAUX_H_

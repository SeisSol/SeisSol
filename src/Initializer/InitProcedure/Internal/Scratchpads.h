// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_SCRATCHPADS_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_SCRATCHPADS_H_

#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer::internal {

/**
 * Derives sizes of scratch memory required during computations of Wave Propagation solver
 **/
void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage);

/**
 * Derives sizes of scratch memory required during computations of Dynamic Rupture solver
 **/
void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_SCRATCHPADS_H_

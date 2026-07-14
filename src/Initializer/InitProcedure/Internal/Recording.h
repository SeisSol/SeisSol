// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_RECORDING_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_RECORDING_H_

#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"

namespace seissol::initializer::internal {

void setupRecorders(LTS::Storage& ltsStorage,
                    DynamicRupture::Storage& drStorage,
                    bool usePlasticity);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_RECORDING_H_

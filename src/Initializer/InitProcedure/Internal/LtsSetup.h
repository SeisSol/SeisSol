// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_

#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Descriptor/LTS.h"
namespace seissol::initializer::internal {
void deriveLtsSetups(const MeshLayout& layout, LTS::Storage& storage);
} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_

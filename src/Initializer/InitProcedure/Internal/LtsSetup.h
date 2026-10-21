// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_

#include "Common/Constants.h"
#include "Initializer/CellLocalInformation.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/TimeStepping/Halo.h"
#include "Memory/Descriptor/LTS.h"

#include <array>
#include <cstdint>

namespace seissol::initializer::internal {

/**
 * Derives the storage requirements of a single cell. Exposed for testing; see the definition for
 * the encoding of the result.
 */
LtsSetup getLtsSetup(const CellLocalInformation& ownPrimary,
                     const SecondaryCellLocalInformation& ownSecondary,
                     const std::array<std::uint64_t, Cell::NumFaces>& neighborClusters);

void deriveLtsSetups(const MeshLayout& layout, LTS::Storage& storage);
} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_LTSSETUP_H_

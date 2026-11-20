// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_CORE_PARALLEL_PARALLELFS_H_
#define SEISSOL_SRC_CORE_PARALLEL_PARALLELFS_H_

#include <filesystem>

namespace seissol {

/**
    Creates a directory, with its parents, for all ranks.

    Collective. Rank 0 creates the directory, then every rank checks once that it can see it.
    Aborts if the directory cannot be created, or if some rank cannot see it.
 */
void createDirectoryGlobally(const std::filesystem::path& path);

} // namespace seissol
#endif // SEISSOL_SRC_CORE_PARALLEL_PARALLELFS_H_

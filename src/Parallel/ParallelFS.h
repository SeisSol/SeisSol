// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_PARALLELFS_H_
#define SEISSOL_SRC_PARALLEL_PARALLELFS_H_

#include "Common/Filesystem.h"

namespace seissol {

/**
    Checks if a file or directory is seen by all ranks in the global communicator.
 */
bool entryExistsGlobally(const filesystem::directory_entry& entry);

/**
    Wait until the path is visible for all ranks.
 */
void waitExistence(const filesystem::path& path);

/**
    Moves a file to a new place, and waits until all ranks see it.
 */
void moveEntryGlobally(const filesystem::path& entry, const filesystem::path& newentry);

/**
    Creates a directory and waits until all ranks see it.
 */
void createDirectoryGlobally(const filesystem::path& entry);

/**
    Turns a file into a backup file if it exists.
 */
void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp);

} // namespace seissol
#endif // SEISSOL_SRC_PARALLEL_PARALLELFS_H_

// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_FILESYSTEM_H_
#define SEISSOL_SRC_COMMON_FILESYSTEM_H_

#include <optional>
#include <string>

#ifdef EXPERIMENTAL_FS
#include <experimental/filesystem> // IWYU pragma: export
#include <sys/stat.h>
#include <sys/types.h>

namespace seissol {
namespace filesystem = std::experimental::filesystem;
} // namespace seissol

#else

#include <filesystem> // IWYU pragma: export
namespace seissol {
namespace filesystem = std::filesystem;
} // namespace seissol

#endif // EXPERIMENTAL_FS

namespace seissol {
auto directoryExists(const seissol::filesystem::directory_entry& entry) -> bool;

/**
    Renames an existing file to a timestamped backup name.

    Purely local: it inspects and renames only files of the calling rank, and does not
    communicate. Not for files which several ranks share.
 */
void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp = {});
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_FILESYSTEM_H_

// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_FILESYSTEM_H_
#define SEISSOL_SRC_COMMON_FILESYSTEM_H_

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

namespace seissol {

/// strftime format of the time stamps in the names of backup files
constexpr auto BackupTimeStampFormat = "%Y-%m-%d_%H-%M-%S";

auto directoryExists(const std::filesystem::directory_entry& entry) -> bool;

/**
    Renames an existing file to a timestamped backup name.

    Purely local: it inspects and renames only files of the calling rank, and does not
    communicate. Not for files which several ranks share.
 */
void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp = {});

/**
 * @brief Replaces the file @p path by one holding @p content, in a single step.
 *
 * The content goes to a file next to it first, which is then renamed to @p path; a reader sees
 * either the old file or the new one, never one that is half written. For files that are
 * rewritten while a run goes on and that other programs may read at any time.
 *
 * @return Whether the file could be written; it is left as it was otherwise.
 */
bool writeFileAtomically(const std::string& path, std::string_view content);
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_FILESYSTEM_H_

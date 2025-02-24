// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
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

void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp = {});
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_FILESYSTEM_H_

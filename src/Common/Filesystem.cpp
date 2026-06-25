// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Filesystem.h"

#include <ctime>
#include <optional>
#include <sstream>
#include <string>
#include <utils/timeutils.h>

namespace seissol {
auto directoryExists(const std::filesystem::directory_entry& entry) -> bool {
  return entry.exists();
}

void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp) {
  std::stringstream fullName;
  fullName << fileName << '.' << fileExtension;
  const auto path = std::filesystem::path(fullName.str());
  const auto entry = std::filesystem::directory_entry(path);

  if (seissol::directoryExists(entry)) {
    const auto actualTimeStamp = timeStamp.value_or(
        []() { return utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr)); }());
    std::stringstream backupFileName;
    backupFileName << fileName << ".bak_" << actualTimeStamp << '.' << fileExtension;
    const auto copyPath = std::filesystem::path(backupFileName.str());
    std::filesystem::rename(path, copyPath);
  }
}
} // namespace seissol

// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Filesystem.h"

#include "utils/timeutils.h"
#include <ctime>
#include <optional>
#include <sstream>
#include <string>

namespace seissol {
#ifdef EXPERIMENTAL_FS
auto directoryExists(const seissol::filesystem::directory_entry& entry) -> bool {
  auto pathName = entry.path().string();
  struct stat info;
  return stat(pathName.c_str(), &info) == 0;
}
#else
auto directoryExists(const seissol::filesystem::directory_entry& entry) -> bool {
  return entry.exists();
}
#endif // EXPERIMENTAL_FS

void generateBackupFileIfNecessary(const std::string& fileName,
                                   const std::string& fileExtension,
                                   const std::optional<std::string>& timeStamp) {
  std::stringstream fullName;
  fullName << fileName << '.' << fileExtension;
  const auto path = seissol::filesystem::path(fullName.str());
  const auto entry = seissol::filesystem::directory_entry(path);

  if (seissol::directoryExists(entry)) {
    const auto actualTimeStap = timeStamp.value_or(
        []() { return utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr)); }());
    std::stringstream backupFileName;
    backupFileName << fileName << ".bak_" << timeStamp.value() << '.' << fileExtension;
    const auto copyPath = seissol::filesystem::path(backupFileName.str());
    seissol::filesystem::rename(path, copyPath);
  }
}
} // namespace seissol

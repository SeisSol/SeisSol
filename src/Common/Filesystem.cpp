// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Filesystem.h"

#include <ctime>
#include <fstream>
#include <ios>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <utils/timeutils.h>

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
    const auto actualTimeStamp = timeStamp.value_or(
        []() { return utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(nullptr)); }());
    std::stringstream backupFileName;
    backupFileName << fileName << ".bak_" << actualTimeStamp << '.' << fileExtension;
    const auto copyPath = seissol::filesystem::path(backupFileName.str());
    seissol::filesystem::rename(path, copyPath);
  }
}

bool writeFileAtomically(const std::string& path, std::string_view content) {
  // next to the target, so that the rename stays within one file system
  const auto temporary = path + ".partial";
  {
    std::ofstream stream(temporary, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!stream.good()) {
      return false;
    }
    stream.write(content.data(), static_cast<std::streamsize>(content.size()));
    stream.flush();
    if (!stream.good()) {
      stream.close();
      std::error_code ignored;
      seissol::filesystem::remove(temporary, ignored);
      return false;
    }
  }
  std::error_code error;
  seissol::filesystem::rename(temporary, path, error);
  if (error) {
    std::error_code ignored;
    seissol::filesystem::remove(temporary, ignored);
    return false;
  }
  return true;
}
} // namespace seissol

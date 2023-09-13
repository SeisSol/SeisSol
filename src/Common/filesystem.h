#ifndef SEISSOL_FILESYSTEM_H
#define SEISSOL_FILESYSTEM_H

#include <string>
#include <optional>
#include "utils/timeutils.h"

#ifdef EXPERIMENTAL_FS
#include <experimental/filesystem>
#include <sys/types.h>
#include <sys/stat.h>

namespace seissol {
namespace filesystem = std::experimental::filesystem;

inline bool directoryExists(seissol::filesystem::directory_entry entry) {
  auto pathName = entry.path().string();
  struct stat info;
  return stat(pathName.c_str(), &info) == 0;
}
} // namespace seissol

#else

#include <filesystem>
namespace seissol {
namespace filesystem = std::filesystem;

inline bool directoryExists(seissol::filesystem::directory_entry entry) {
  return entry.exists();
}
} // namespace seissol

#endif // EXPERIMENTAL_FS


namespace seissol {
inline void generateBackupFileIfNecessary(std::string fileName,
                                          std::string fileExtension,
                                          std::optional<std::string> timeStamp = {}) {
  std::stringstream fullName;
  fullName << fileName << '.' << fileExtension;
  seissol::filesystem::path path(fullName.str());
  seissol::filesystem::directory_entry entry(path);

  if (seissol::directoryExists(entry)) {
    if (!timeStamp.has_value()) {
      auto stamp = utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", time(0L));
      timeStamp = std::optional<std::string>(stamp);
    }
    std::stringstream backupFileName;
    backupFileName << fileName << ".bak_" << timeStamp.value() << '.' << fileExtension;
    seissol::filesystem::path copyPath(backupFileName.str());
    seissol::filesystem::rename(path, copyPath);
  }
}
} // namespace seissol

#endif // SEISSOL_FILESYSTEM_H


#ifndef SEISSOL_FILESYSTEM_H
#define SEISSOL_FILESYSTEM_H

#ifdef EXPERIMENTAL_FS

#include <experimental/filesystem>
#include <sys/types.h>
#include <sys/stat.h>

namespace seissol {
namespace filesystem = std::experimental::filesystem;

inline bool directoryExists(seissol::filesystem::directory_entry entry) {
  auto pathName = entry.path().string();
  struct stat info;
  return stat(pathName.c_str(), &info) != 0;
}
} // namespace seissol

#else

#include <filesystem>
namespace seissol {
namespace filesystem = std::filesystem;

inline bool directoryExists(seissol::filesystem::directory_entry entry) { return entry.exists(); }
} // namespace seissol

#endif // EXPERIMENTAL_FS

#endif // SEISSOL_FILESYSTEM_H

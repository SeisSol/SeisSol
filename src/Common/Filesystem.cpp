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
} // namespace seissol

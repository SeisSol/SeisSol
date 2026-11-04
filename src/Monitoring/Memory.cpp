// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Memory.h"

#include <cstddef>
#include <optional>

#ifndef _WIN32
#include <sys/resource.h>
#endif

namespace seissol::monitoring {

std::optional<std::size_t> peakHostMemory() {
#ifdef _WIN32
  return {};
#else
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return {};
  }
  // glibc puts every field of rusage in a union with a word of the kernel interface
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
  const auto maxResident = static_cast<std::size_t>(usage.ru_maxrss);
#ifdef __APPLE__
  // in bytes on macOS
  return maxResident;
#else
  // in KiB on Linux and the BSDs
  return maxResident * 1024;
#endif
#endif
}

} // namespace seissol::monitoring

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Environment.h"

// setenv and unsetenv are POSIX and only declared by the C header; Windows has _putenv_s instead
#include <stdlib.h> // NOLINT(modernize-deprecated-headers)

namespace seissol {

void setEnvironmentVariable(const char* name, const char* value) noexcept {
#ifdef _WIN32
  _putenv_s(name, value);
#else
  // NOLINTNEXTLINE(concurrency-mt-unsafe)
  setenv(name, value, 1);
#endif
}

void unsetEnvironmentVariable(const char* name) noexcept {
#ifdef _WIN32
  // (an empty value removes the variable)
  _putenv_s(name, "");
#else
  // NOLINTNEXTLINE(concurrency-mt-unsafe)
  unsetenv(name);
#endif
}

} // namespace seissol

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Library.h"

#include <cstdlib>
#include <dlfcn.h>
#include <optional>
#include <string>
#include <vector>

namespace seissol::monitoring::energy {

std::optional<Library> Library::open(const std::vector<std::string>& names) {
  for (const auto& name : names) {
    // NOLINTNEXTLINE(hicpp-signed-bitwise)
    void* handle = dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (handle != nullptr) {
      return Library(handle, name);
    }
  }
  return {};
}

void* Library::symbol(const char* name) const { return dlsym(handle_, name); }

std::vector<std::string> rocmLibraryNames(const std::string& name) {
  std::vector<std::string> names{name};
  // NOLINTNEXTLINE(concurrency-mt-unsafe)
  if (const char* rocm = std::getenv("ROCM_PATH"); rocm != nullptr && *rocm != '\0') {
    names.push_back(std::string(rocm) + "/lib/" + name);
  }
  names.push_back("/opt/rocm/lib/" + name);
  return names;
}

} // namespace seissol::monitoring::energy

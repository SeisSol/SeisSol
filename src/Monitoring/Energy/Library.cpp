// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Library.h"

#include <cstdlib>
#include <optional>
#include <string>
#include <vector>

#ifndef _WIN32
#include <dlfcn.h>
#endif

namespace seissol::monitoring::energy {

std::optional<Library> Library::open(const std::vector<std::string>& names) {
#ifdef _WIN32
  // the names are those of Linux; the nodes that SeisSol measures the energy of run it
  static_cast<void>(names);
#else
  for (const auto& name : names) {
    // NOLINTNEXTLINE(hicpp-signed-bitwise)
    void* handle = dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (handle != nullptr) {
      return Library(handle, name);
    }
  }
#endif
  return {};
}

void* Library::symbol(const char* name) const {
#ifdef _WIN32
  static_cast<void>(name);
  return nullptr;
#else
  return dlsym(handle_, name);
#endif
}

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

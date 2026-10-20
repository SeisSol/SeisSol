// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostData.h"

#include <cassert>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::initializer {

HostData HostData::everything() {
  HostData data;
  data.everything_ = true;
  return data;
}

HostData& HostData::merge(const HostData& other) {
  everything_ = everything_ || other.everything_;
  for (const auto& [variable, copies] : other.variables_) {
    variables_.try_emplace(variable, copies);
  }
  return *this;
}

std::vector<HostCopy> HostData::copies() const {
  assert(!everything_);
  std::vector<HostCopy> copies;
  for (const auto& [key, addCopies] : variables_) {
    const auto& [storage, variable] = key;
    // a single variable of a storage whose variables all get copied anyway
    if (variable != AllVariables && variables_.count(std::make_pair(storage, AllVariables)) > 0) {
      continue;
    }
    addCopies(copies);
  }
  return copies;
}

void HostData::copyToHost([[maybe_unused]] void* stream) const {
#ifdef ACL_DEVICE
  auto& api = device::DeviceInstance::instance().api();
  for (const auto& copy : copies()) {
    api.copyFromAsync(copy.host, copy.device, copy.size, stream);
  }
#endif
}

} // namespace seissol::initializer

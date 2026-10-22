// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostData.h"

#include "Memory/Tree/Layer.h"

#include <cassert>
#include <utility>
#include <vector>

namespace seissol::initializer {

HostData HostData::everything() {
  HostData data;
  data.everything_ = true;
  return data;
}

HostData& HostData::merge(const HostData& other) {
  everything_ = everything_ || other.everything_;
  for (const auto& [variable, containers] : other.variables_) {
    variables_.try_emplace(variable, containers);
  }
  return *this;
}

std::vector<const DualMemoryContainer*> HostData::containers() const {
  assert(!everything_);
  std::vector<const DualMemoryContainer*> containers;
  for (const auto& [key, addContainers] : variables_) {
    const auto& [storage, variable] = key;
    // a single variable of a storage whose variables all get copied anyway
    if (variable != AllVariables && variables_.count(std::make_pair(storage, AllVariables)) > 0) {
      continue;
    }
    addContainers(containers);
  }
  return containers;
}

void HostData::copyToHost(void* stream) const {
  for (const auto* container : containers()) {
    container->synchronizeTo(AllocationPlace::Host, stream);
  }
}

} // namespace seissol::initializer

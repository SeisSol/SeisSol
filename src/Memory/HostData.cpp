// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HostData.h"

#include <cassert>
#include <utility>

namespace seissol::initializer {

HostData HostData::everything() {
  HostData data;
  data.everything_ = true;
  return data;
}

HostData& HostData::merge(const HostData& other) {
  everything_ = everything_ || other.everything_;
  for (const auto& [variable, copy] : other.copies_) {
    copies_.try_emplace(variable, copy);
  }
  return *this;
}

void HostData::copyToHost(void* stream) const {
  assert(!everything_);
  for (const auto& [key, copy] : copies_) {
    const auto& [storage, variable] = key;
    // a single variable of a storage whose variables all get copied anyway
    if (variable != AllVariables && copies_.count(std::make_pair(storage, AllVariables)) > 0) {
      continue;
    }
    copy(stream);
  }
}

} // namespace seissol::initializer

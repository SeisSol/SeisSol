// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "RunConfigs.h"

#include "Core/Common/ConfigRegistry.h"

#include <algorithm>
#include <vector>

namespace seissol::initializer::parameters {

ConfigId RunConfigs::configOfGroup(int group) const {
  const auto found = groupConfigs.find(group);
  return found == groupConfigs.end() ? config : found->second;
}

std::vector<ConfigId> RunConfigs::configs() const {
  std::vector<ConfigId> ofGroups;
  ofGroups.reserve(groupConfigs.size());
  for (const auto& [group, groupConfig] : groupConfigs) {
    ofGroups.push_back(groupConfig);
  }
  std::ranges::sort(ofGroups);
  std::vector<ConfigId> result{config};
  for (const auto groupConfig : ofGroups) {
    if (std::ranges::find(result, groupConfig) == result.end()) {
      result.push_back(groupConfig);
    }
  }
  return result;
}

} // namespace seissol::initializer::parameters

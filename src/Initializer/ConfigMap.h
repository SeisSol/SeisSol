// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_CONFIGMAP_H_
#define SEISSOL_SRC_INITIALIZER_CONFIGMAP_H_

#include "Config.h"

#include <Common/ConfigHelper.h>
#include <Common/Iterator.h>
#include <unordered_map>
#include <utils/env.h>

namespace seissol {

class ConfigMap {
  public:
  ConfigMap(const std::unordered_map<int, std::string>& input, utils::Env& env) {
    bool hasDefault = false;
    const auto defaultCfgString = env.getOptional<std::string>("DEFAULT_CONFIG");
    if (defaultCfgString.has_value()) {
      for (auto [i, str] : common::enumerate(ConfigString)) {
        if (str == defaultCfgString.value()) {
          defaultCfg = ConfigVariantList[i];
          hasDefault = true;
          break;
        }
      }
      if (!hasDefault) {
        logWarning() << "Default config" << defaultCfgString.value()
                     << "given, but not found in the current executable. It only has"
                     << ConfigString;
      }
    }

    if (!hasDefault) {
      defaultCfg = ConfigVariantList[0];
    }

    for (const auto& [group, configString] : input) {
      for (auto [i, str] : common::enumerate(ConfigString)) {
        if (str == defaultCfgString.value()) {
          map[group] = ConfigVariantList[i];
          break;
        }
      }
    }
  }

  [[nodiscard]] ConfigVariant toConfig(int group) const {
    auto mapSearch = map.find(group);
    if (mapSearch == map.end()) {
      return defaultCfg;
    } else {
      return mapSearch->second;
    }
  }

  private:
  std::unordered_map<int, ConfigVariant> map;
  ConfigVariant defaultCfg;
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_CONFIGMAP_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "FaceMap.h"

#include "Initializer/BasicTypedefs.h"

#include <utils/logger.h>
#include <utils/stringutils.h>
#include <yaml-cpp/yaml.h>

namespace seissol {

FaceMap parseFaceMap(YAML::Node node) {
  FaceMap map;

  const static std::unordered_map<std::string, FaceType> stringToNameMap = {
      {"regular", FaceType::Regular},
      {"freeSurface", FaceType::FreeSurface},
      {"freeSurfaceGravity", FaceType::FreeSurfaceGravity},
      {"dynamicRupture", FaceType::DynamicRupture},
      {"dirichlet", FaceType::Dirichlet},
      {"outflow", FaceType::Outflow},
      {"analytic", FaceType::Analytical},
  };

  for (const auto& entry : node) {
    const auto stringType = entry.first.as<std::string>();
    const auto typeFind = stringToNameMap.find(stringType);
    if (typeFind == stringToNameMap.end()) {
      logError() << "";
    }
    const auto type = typeFind->second;

    const auto processEntry = [&](const auto& listEntry) {
      auto strParts = utils::StringUtils::split(listEntry, ',');
      for (auto& part : strParts) {
        std::optional<uint32_t> lower;
        std::optional<uint32_t> upper;

        utils::StringUtils::trim(part);
        auto subParts = utils::StringUtils::split(part, '-');
        for (auto& subPart : subParts) {
          utils::StringUtils::trim(subPart);
        }
        if (subParts.size() == 1) {
          lower = std::stoi(subParts[0]);
          upper = lower;
        } else if (subParts.size() == 2) {
          if (subParts[0].empty()) {
            lower = {};
          } else {
            lower = std::stoi(subParts[0]);
          }
          if (subParts[1].empty()) {
            upper = {};
          } else {
            upper = std::stoi(subParts[0]);
          }
        } else {
          logError() << "Parsing error.";
        }

        map.addRange(lower, upper, type);
      }
    };

    if (entry.second.IsSequence()) {
      for (const auto& listEntry : entry.second) {
        processEntry(listEntry.as<std::string>());
      }
    } else {
      processEntry(entry.second.as<std::string>());
    }
  }

  return map;
}

FaceMap defaultFaceMap() {
  FaceMap map;

  map.addEntry(0, FaceType::Regular);
  map.addEntry(6, FaceType::Regular);

  map.addEntry(1, FaceType::FreeSurface);
  map.addEntry(2, FaceType::FreeSurfaceGravity);

  map.addEntry(3, FaceType::DynamicRupture);
  map.addRange(65, {}, FaceType::DynamicRupture);

  map.addEntry(4, FaceType::Dirichlet);

  map.addEntry(5, FaceType::Outflow);

  map.addEntry(7, FaceType::Analytical);

  return map;
}

} // namespace seissol

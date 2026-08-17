// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "FaceMap.h"

#include "Initializer/BasicTypedefs.h"

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <yaml-cpp/yaml.h>

namespace seissol {

FaceMap parseFaceMap(const YAML::Node& node) {
  FaceMap map;

  const static std::unordered_map<std::string, FaceType> StringToNameMap = {
      {"regular", FaceType::Regular},
      {"freeSurface", FaceType::FreeSurface},
      {"freeSurfaceGravity", FaceType::FreeSurfaceGravity},
      {"dynamicRupture", FaceType::DynamicRupture},
      {"dirichlet", FaceType::Dirichlet},
      {"outflow", FaceType::Outflow},
      {"analytical", FaceType::Analytical},
  };

  for (const auto& entry : node) {
    const auto stringType = entry.first.as<std::string>();
    const auto typeFind = StringToNameMap.find(stringType);
    if (typeFind == StringToNameMap.end()) {
      logError() << "";
    }
    const auto type = typeFind->second;

    const auto processEntry = [&](const auto& listEntry) {
      auto strParts = utils::StringUtils::split(listEntry, ',');
      if (utils::StringUtils::endsWith(listEntry, ",")) {
        // FIXME: the StringUtils::split does not consider TRAILING '-'s in this case.
        strParts.emplace_back();
      }

      for (auto& part : strParts) {
        std::optional<uint32_t> lower;
        std::optional<uint32_t> upper;

        utils::StringUtils::trim(part);
        auto subParts = utils::StringUtils::split(part, '-');

        if (utils::StringUtils::endsWith(part, "-")) {
          // FIXME: the StringUtils::split does not consider TRAILING '-'s in this case.
          subParts.emplace_back();
        }

        for (auto& subPart : subParts) {
          utils::StringUtils::trim(subPart);
        }
        if (subParts.size() == 1) {
          lower = std::stoi(subParts[0]);
          upper = lower;
        } else if (subParts.size() == 2) {
          if (subParts[0] == "") {
            lower = {};
          } else {
            lower = std::stoi(subParts[0]);
          }
          if (subParts[1] == "") {
            upper = {};
          } else {
            upper = std::stoi(subParts[1]);
          }
        } else {
          logError() << "Parsing error while parsing face map. At:" << part << "i.e." << subParts;
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

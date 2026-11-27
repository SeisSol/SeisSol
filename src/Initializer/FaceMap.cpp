// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "FaceMap.h"

#include "Common/SegmentMap.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <optional>
#include <sstream>
#include <string>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol {

namespace {

// the ranges of a comma-separated list such as "1, 3-5, 65-": a number, or a range whose bounds
// are inclusive and may be open
template <typename F>
void forEachRange(const std::string& list, F&& add) {
  auto strParts = utils::StringUtils::split(list, ',');
  if (utils::StringUtils::endsWith(list, ",")) {
    // FIXME: the StringUtils::split does not consider TRAILING '-'s in this case.
    strParts.emplace_back();
  }

  for (auto& part : strParts) {
    std::optional<std::int64_t> lower;
    std::optional<std::int64_t> upper;

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

    add(lower, upper);
  }
}

// the ranges of a YAML value: such a list, or a sequence of them
template <typename F>
void forEachRangeOf(const YAML::Node& value, F&& add) {
  if (value.IsSequence()) {
    for (const auto& listEntry : value) {
      forEachRange(listEntry.as<std::string>(), add);
    }
  } else {
    forEachRange(value.as<std::string>(), add);
  }
}

// a range as the keys of a segment map take it
template <typename KeyT>
std::optional<KeyT> boundOf(const std::optional<std::int64_t>& bound) {
  if (!bound.has_value()) {
    return {};
  }
  return static_cast<KeyT>(bound.value());
}

// the entries of a face type: a list of tags, or a sequence of such lists and of maps, which
// name their tags under "tags"
std::vector<YAML::Node> entriesOf(const YAML::Node& value) {
  std::vector<YAML::Node> entries;
  if (value.IsSequence()) {
    for (const auto& listEntry : value) {
      entries.push_back(listEntry);
    }
  } else {
    entries.push_back(value);
  }
  return entries;
}

FaceType faceTypeOf(const std::string& name) {
  const auto* const typeFind = std::ranges::find_if(
      FaceTypes, [&](FaceType candidate) { return faceTypeName(candidate) == name; });
  if (typeFind == FaceTypes.end()) {
    std::stringstream known;
    for (std::size_t i = 0; i < FaceTypes.size(); ++i) {
      known << (i == 0 ? "" : ", ") << faceTypeName(FaceTypes[i]);
    }
    logError() << "Unknown face type" << name << "in the face map. Known types are:" << known.str();
  }
  return *typeFind;
}

// the text of a list of tags or groups, for the log
std::string listText(const YAML::Node& value) {
  if (!value.IsSequence()) {
    return value.as<std::string>();
  }
  std::string text;
  for (const auto& listEntry : value) {
    text += (text.empty() ? "" : ", ") + listEntry.as<std::string>();
  }
  return text;
}

} // namespace

FaceMap parseFaceMap(const YAML::Node& node) {
  FaceMap map;

  for (const auto& entry : node) {
    const auto type = faceTypeOf(entry.first.as<std::string>());
    for (const auto& listEntry : entriesOf(entry.second)) {
      const auto tags = listEntry.IsMap() ? listEntry["tags"] : listEntry;
      if (tags) {
        forEachRangeOf(tags, [&](const auto& lower, const auto& upper) {
          map.addRange(boundOf<std::uint32_t>(lower), boundOf<std::uint32_t>(upper), type);
        });
      }
    }
  }

  return map;
}

std::vector<initializer::parameters::FaceFluxRule> parseFaceFluxRules(const YAML::Node& node) {
  using initializer::parameters::FaceFluxRule;
  std::vector<FaceFluxRule> rules;

  for (const auto& entry : node) {
    const auto typeName = entry.first.as<std::string>();
    const auto type = faceTypeOf(typeName);
    std::size_t index = 0;
    for (const auto& listEntry : entriesOf(entry.second)) {
      ++index;
      if (!listEntry.IsMap()) {
        continue;
      }
      if (type != FaceType::Regular) {
        for (const auto& field : listEntry) {
          const auto key = field.first.as<std::string>();
          if (key != "tags") {
            logError() << "The face map gives" << key.c_str() << "for" << typeName.c_str()
                       << "faces; only regular faces take more than their tags.";
          }
        }
        continue;
      }

      FaceFluxRule rule;
      rule.name = "the entry " + std::to_string(index) + " of regular in the face map";
      if (listEntry["tags"]) {
        rule.name += ", tags " + listText(listEntry["tags"]);
        SegmentMap<std::uint32_t, bool> tags;
        forEachRangeOf(listEntry["tags"], [&](const auto& lower, const auto& upper) {
          tags.addRange(boundOf<std::uint32_t>(lower), boundOf<std::uint32_t>(upper), true);
        });
        rule.tags = tags;
      }
      if (listEntry["groups"]) {
        rule.name += ", groups " + listText(listEntry["groups"]);
        SegmentMap<int, bool> groups;
        forEachRangeOf(listEntry["groups"], [&](const auto& lower, const auto& upper) {
          groups.addRange(boundOf<int>(lower), boundOf<int>(upper), true);
        });
        rule.groups = groups;
      }

      initializer::parameters::ParameterReader reader(listEntry, "", false);
      reader.markUnused({"tags", "groups"});
      const auto setsFlux = [&](const std::string& suffix) {
        return reader.hasField("numflux" + suffix) || reader.hasField("fluxgpart" + suffix) ||
               reader.hasField("fluxrpart" + suffix) || reader.hasField("fluxdamping" + suffix);
      };
      if (setsFlux("")) {
        rule.flux = initializer::parameters::readFluxSetting(&reader, "");
      }
      if (setsFlux("nearfault")) {
        rule.fluxNearFault = initializer::parameters::readFluxSetting(&reader, "nearfault");
      }
      reader.warnUnknown(rule.name);

      if (!rule.flux.has_value() && !rule.fluxNearFault.has_value()) {
        if (rule.groups.has_value()) {
          logError() << "The face map selects faces by their mesh groups in" << rule.name.c_str()
                     << "but gives them no flux.";
        }
        // tags alone, written as a map
        continue;
      }
      if (!rule.tags.has_value() && !rule.groups.has_value()) {
        logError() << "The face map gives a flux in" << rule.name.c_str()
                   << "but selects no faces: it takes tags, groups, or both. The flux of all "
                      "regular faces is the one of the parameter file.";
      }
      rules.push_back(rule);
    }
  }

  return rules;
}

std::vector<initializer::parameters::FaceFluxRule>
    readFaceFluxRules(const std::optional<std::string>& path) {
  if (!path.has_value()) {
    return {};
  }
  YAML::Node node;
  try {
    node = YAML::LoadFile(path.value());
  } catch (const std::exception& error) {
    logError() << "Could not read the face map file" << path.value() << ":"
               << std::string(error.what());
  }
  return parseFaceFluxRules(node);
}

FaceMap readFaceMap(const std::string& path) {
  YAML::Node node;
  try {
    node = YAML::LoadFile(path);
  } catch (const std::exception& error) {
    logError() << "Could not read the face map file" << path << ":" << std::string(error.what());
  }
  return parseFaceMap(node);
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

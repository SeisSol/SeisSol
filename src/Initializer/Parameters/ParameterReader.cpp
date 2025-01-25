// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ParameterReader.h"
#include <memory>
#include <optional>
#include <string>
#include <utils/logger.h>
#include <vector>

#include "Common/Filesystem.h"

namespace seissol::initializer::parameters {

ParameterReader::ParameterReader(const YAML::Node& node, const std::string& rootPath, bool empty)
    : node(node), rootPath(rootPath), empty(empty) {}

std::optional<std::string> ParameterReader::readPath(const std::string& field) {
  const auto fileName = read<std::string>(field);
  if (fileName.has_value()) {
    const auto lastPath = filesystem::path(rootPath);
    auto nextPath = filesystem::path(fileName.value());
    const auto loadFileName = [&]() {
      if (nextPath.is_relative()) {
        // remove file
        return lastPath.parent_path() / nextPath;
      } else {
        return nextPath;
      }
    }();
    // try to get the full file name (if that works)
    if (filesystem::exists(loadFileName)) {
      return filesystem::canonical(loadFileName);
    } else if (filesystem::exists(nextPath)) {
      return filesystem::canonical(nextPath);
    } else {
      // otherwise, just return the string (TODO: or should be fail here then?)
      return nextPath;
    }
  } else {
    return {};
  }
}

std::string ParameterReader::readPathOrFail(const std::string& field,
                                            const std::string& failMessage) {
  const auto path = readPath(field);
  if (path.has_value()) {
    return path.value();
  } else {
    logError() << "The field" << field
               << "was not found, but it is required (details: " << failMessage.c_str() << ")";
    return "";
  }
}

void ParameterReader::warnDeprecatedSingle(const std::string& field) {
  if (hasField(field)) {
    visited.emplace(field);
    logInfo() << "The field" << field
              << "is no longer in use. You may safely remove it from your parameters file.";
  }
}

void ParameterReader::warnDeprecated(const std::vector<std::string>& fields) {
  for (const auto& field : fields) {
    warnDeprecatedSingle(field);
  }
}

void ParameterReader::warnUnknown(const std::string& prefix) const {
  for (const auto& subnodes : node) {
    auto field = subnodes.first.as<std::string>();
    if (visited.find(field) == visited.end()) {
      logWarning() << "The field" << field << "in" << prefix
                   << "was given in the parameter file, but is unknown to SeisSol.";
    }
  }
  for (const auto& pair : subreaders) {
    pair.second->warnUnknown(pair.first);
  }
}

void ParameterReader::markUnused(const std::vector<std::string>& fields) {
  for (const auto& field : fields) {
    logDebug() << "The field" << field << "is ignored (if it is found).";
    visited.emplace(field);
  }
}

ParameterReader* ParameterReader::readSubNode(const std::string& subnodeName) {
  visited.emplace(subnodeName);
  logDebug() << "Entering section" << subnodeName;
  if (subreaders.find(subnodeName) == subreaders.end()) {
    bool empty = false;
    if (hasField(subnodeName)) {
      empty = false;

    } else {
      logDebug() << "Section" << subnodeName
                 << "not found in the given parameter file. Using an empty reader->";
      empty = true;
    }
    subreaders.emplace(
        subnodeName,
        std::make_shared<ParameterReader>(ParameterReader(node[subnodeName], rootPath, empty)));
  }
  return subreaders.at(subnodeName).get();
}

bool ParameterReader::hasField(const std::string& field) { return !empty && node[field]; }

} // namespace seissol::initializer::parameters

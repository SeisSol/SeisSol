// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ParameterReader.h"

#include "Common/Filesystem.h"

#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::parameters {

ParameterReader::ParameterReader(const YAML::Node& node, const std::string& rootPath, bool empty)
    : node_(node), rootPath_(rootPath), empty_(empty) {}

std::optional<std::string> ParameterReader::readPath(const std::string& field) {
  const auto fileName = read<std::string>(field);
  if (!fileName.has_value()) {
    return {};
  }
  auto resolved = [&]() -> std::string {
    if (fileName->empty()) {
      // no file; relative to the directory of the parameter file, it would name that directory
      return fileName.value();
    }
    const auto lastPath = std::filesystem::path(rootPath_);
    auto nextPath = std::filesystem::path(fileName.value());
    const auto loadFileName = [&] {
      if (nextPath.is_relative()) {
        // remove file
        return lastPath.parent_path() / nextPath;
      } else {
        return nextPath;
      }
    }();
    // try to get the full file name (if that works)
    if (std::filesystem::exists(loadFileName)) {
      return std::filesystem::canonical(loadFileName);
    } else if (std::filesystem::exists(nextPath)) {
      return std::filesystem::canonical(nextPath);
    } else {
      // otherwise, just return the string (TODO: or should be fail here then?)
      return nextPath;
    }
  }();
  paths_.emplace_back(field, resolved);
  return resolved;
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
    visited_.emplace(field);
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
  for (const auto& subnodes : node_) {
    const auto field = subnodes.first.as<std::string>();
    if (!visited_.contains(field)) {
      logWarning() << "The field" << field << "in" << prefix
                   << "was given in the parameter file, but is unknown to SeisSol.";
    }
  }
  for (const auto& pair : subreaders_) {
    pair.second->warnUnknown(pair.first);
  }
}

void ParameterReader::markUnused(const std::vector<std::string>& fields) {
  for (const auto& field : fields) {
    logDebug() << "The field" << field << "is ignored (if it is found).";
    visited_.emplace(field);
  }
}

ParameterReader* ParameterReader::readSubNode(const std::string& subnodeName) {
  visited_.emplace(subnodeName);
  logDebug() << "Entering section" << subnodeName;
  if (!subreaders_.contains(subnodeName)) {
    bool empty = false;
    if (hasField(subnodeName)) {
      empty = false;

    } else {
      logDebug() << "Section" << subnodeName
                 << "not found in the given parameter file. Using an empty reader->";
      empty = true;
    }
    subreaders_.emplace(
        subnodeName,
        std::make_shared<ParameterReader>(ParameterReader(node_[subnodeName], rootPath_, empty)));
  }
  return subreaders_.at(subnodeName).get();
}

bool ParameterReader::hasField(const std::string& field) { return !empty_ && node_[field]; }

namespace {
//! @brief The names of the sections in @p subreaders , in alphabetical order.
template <typename Map>
std::vector<std::string> sortedSections(const Map& subreaders) {
  std::vector<std::string> sections;
  sections.reserve(subreaders.size());
  for (const auto& [name, _] : subreaders) {
    sections.push_back(name);
  }
  std::ranges::sort(sections);
  return sections;
}
} // namespace

std::vector<ParameterReader::NamedPath> ParameterReader::paths() const {
  std::vector<NamedPath> result;
  result.reserve(paths_.size());
  for (const auto& [field, path] : paths_) {
    result.push_back(NamedPath{.field = {field}, .path = path});
  }
  for (const auto& section : sortedSections(subreaders_)) {
    for (auto entry : subreaders_.at(section)->paths()) {
      entry.field.insert(entry.field.begin(), section);
      result.push_back(std::move(entry));
    }
  }
  return result;
}

std::vector<std::string> ParameterReader::unknownFields() const {
  std::vector<std::string> result;
  for (const auto& subnodes : node_) {
    const auto field = subnodes.first.as<std::string>();
    if (!visited_.contains(field)) {
      result.push_back(field);
    }
  }
  for (const auto& section : sortedSections(subreaders_)) {
    for (const auto& field : subreaders_.at(section)->unknownFields()) {
      result.emplace_back(section).append(".").append(field);
    }
  }
  return result;
}

bool ParameterReader::visited(const std::vector<std::string>& path) const {
  if (path.empty()) {
    return false;
  }
  if (path.size() == 1) {
    return visited_.contains(path.front());
  }
  const auto subreader = subreaders_.find(path.front());
  if (subreader == subreaders_.end()) {
    return false;
  }
  return subreader->second->visited(std::vector<std::string>(path.begin() + 1, path.end()));
}

} // namespace seissol::initializer::parameters

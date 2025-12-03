// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Hdf5.h"

#include "Data.h"
#include "IO/Datatype/Datatype.h"

#include <cassert>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <utils/stringutils.h>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {
Hdf5Location::Hdf5Location(const std::string& longstring) {
  auto parts = utils::StringUtils::split(longstring, ':');
  assert(parts.size() == 2);
  fileP = parts[0];
  auto groupParts = utils::StringUtils::split(parts[1], '/');
  groupsP = groupParts;
}

Hdf5Location::Hdf5Location(const std::string& file,
                           const std::vector<std::string>& groups,
                           const std::optional<std::string>& dataset)
    : fileP(file), groupsP(groups), datasetP(dataset) {}

Hdf5Location::Hdf5Location(YAML::Node node)
    : fileP(node["file"].as<std::string>()), groupsP(node["group"].as<std::vector<std::string>>()),
      datasetP(node["dataset"] ? node["dataset"].as<std::string>() : std::optional<std::string>()) {
}

std::string Hdf5Location::file() const { return fileP; }
std::vector<std::string> Hdf5Location::groups() const { return groupsP; }
std::optional<std::string> Hdf5Location::dataset() const { return datasetP; }
std::string Hdf5Location::infilePath() const {
  std::string path;
  for (const auto& group : groupsP) {
    path += "/" + group;
  }
  if (datasetP.has_value()) {
    path += "/" + datasetP.value();
  }
  return path;
}

std::optional<Hdf5Location> Hdf5Location::commonLocation(const Hdf5Location& other) const {
  if (other.file() == file()) {
    std::vector<std::string> commonGroups;
    for (std::size_t i = 0; i < groups().size() && i < other.groups().size(); ++i) {
      if (groups()[i] == other.groups()[i]) {
        commonGroups.push_back(groups()[i]);
      } else {
        break;
      }
    }
    return std::make_optional<Hdf5Location>(file(), commonGroups);
  } else {
    return std::optional<Hdf5Location>();
  }
}

YAML::Node Hdf5Location::serialize() {
  YAML::Node node;
  node["file"] = fileP;
  node["group"] = groupsP;
  if (datasetP.has_value()) {
    node["dataset"] = datasetP.value();
  }
  return node;
}

YAML::Node Hdf5AttributeWrite::serialize() {
  YAML::Node node;
  node["name"] = name;
  node["source"] = dataSource->serialize();
  node["location"] = location.serialize();
  node["writer"] = "hdf5";
  node["type"] = "attribute";
  return node;
}

Hdf5AttributeWrite::Hdf5AttributeWrite(const Hdf5Location& location,
                                       const std::string& name,
                                       std::shared_ptr<writer::DataSource> dataSource)
    : location(location), name(name), dataSource(std::move(dataSource)) {}

Hdf5AttributeWrite::Hdf5AttributeWrite(YAML::Node node)
    : location(Hdf5Location(node["location"])), name(node["name"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])) {}

std::vector<std::shared_ptr<DataSource>> Hdf5AttributeWrite::dataSources() { return {dataSource}; }

Hdf5DataWrite::Hdf5DataWrite(const Hdf5Location& location,
                             const std::string& name,
                             std::shared_ptr<writer::DataSource> dataSource,
                             std::shared_ptr<datatype::Datatype> targetType,
                             bool append,
                             int compress)
    : location(location), name(name), dataSource(std::move(dataSource)),
      targetType(std::move(targetType)), append(append), compress(compress) {}

YAML::Node Hdf5DataWrite::serialize() {
  YAML::Node node;
  node["name"] = name;
  node["source"] = dataSource->serialize();
  node["location"] = location.serialize();
  node["targetType"] = targetType->serialize();
  node["writer"] = "hdf5";
  node["type"] = "data";
  node["append"] = append;
  node["compress"] = compress;
  return node;
}

Hdf5DataWrite::Hdf5DataWrite(YAML::Node node)
    : location(Hdf5Location(node["location"])), name(node["name"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])),
      targetType(datatype::Datatype::deserialize(node["targetType"])),
      append(node["append"].as<bool>()), compress(node["compress"].as<int>()) {}

std::vector<std::shared_ptr<DataSource>> Hdf5DataWrite::dataSources() { return {dataSource}; }

YAML::Node Hdf5LinkExternalWrite::serialize() {
  YAML::Node node;
  node["name"] = name;
  node["location"] = location.serialize();
  node["remote"] = remote.serialize();
  node["writer"] = "hdf5";
  node["type"] = "link-external";
  return node;
}

Hdf5LinkExternalWrite::Hdf5LinkExternalWrite(const Hdf5Location& location,
                                             const std::string& name,
                                             const Hdf5Location& remote)
    : location(location), name(name), remote(remote) {}

Hdf5LinkExternalWrite::Hdf5LinkExternalWrite(YAML::Node node)
    : name(node["name"].as<std::string>()), location(Hdf5Location(node["location"])),
      remote(Hdf5Location(node["remote"])) {}

std::vector<std::shared_ptr<DataSource>> Hdf5LinkExternalWrite::dataSources() { return {}; }

} // namespace seissol::io::writer::instructions

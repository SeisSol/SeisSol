// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Hdf5.h"

#include "Data.h"
#include <IO/Datatype/Datatype.h>
#include <cassert>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "utils/stringutils.h"

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
    : name(node["name"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])),
      location(Hdf5Location(node["location"])) {}

std::vector<std::shared_ptr<DataSource>> Hdf5AttributeWrite::dataSources() { return {dataSource}; }

Hdf5DataWrite::Hdf5DataWrite(const Hdf5Location& location,
                             const std::string& name,
                             std::shared_ptr<writer::DataSource> dataSource,
                             std::shared_ptr<datatype::Datatype> targetType,
                             int compress)
    : location(location), name(name), dataSource(std::move(dataSource)),
      targetType(std::move(targetType)), compress(compress) {}

YAML::Node Hdf5DataWrite::serialize() {
  YAML::Node node;
  node["name"] = name;
  node["source"] = dataSource->serialize();
  node["location"] = location.serialize();
  node["targetType"] = targetType->serialize();
  node["writer"] = "hdf5";
  node["type"] = "data";
  node["compress"] = compress;
  return node;
}

Hdf5DataWrite::Hdf5DataWrite(YAML::Node node)
    : name(node["name"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])),
      location(Hdf5Location(node["location"])),
      targetType(datatype::Datatype::deserialize(node["targetType"])),
      compress(node["compress"].as<int>()) {}

std::vector<std::shared_ptr<DataSource>> Hdf5DataWrite::dataSources() { return {dataSource}; }
} // namespace seissol::io::writer::instructions

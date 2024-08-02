#pragma once

#include "Data.hpp"
#include "Instruction.hpp"
#include <memory>
#include <optional>
#include <string>
#include <yaml-cpp/yaml.h>

#include "utils/stringutils.h"

namespace seissol::io::writer::instructions {
class Hdf5Location {
  public:
  Hdf5Location(const std::string& longstring) {
    auto parts = utils::StringUtils::split(longstring, ':');
    assert(parts.size() == 2);
    fileP = parts[0];
    auto groupParts = utils::StringUtils::split(parts[1], '/');
    groupsP = groupParts;
  }

  Hdf5Location(const std::string& file, const std::vector<std::string>& groups)
      : fileP(file), groupsP(groups) {}

  explicit Hdf5Location(YAML::Node node)
      : fileP(node["file"].as<std::string>()),
        groupsP(node["group"].as<std::vector<std::string>>()) {}

  std::string file() const { return fileP; }
  std::vector<std::string> groups() const { return groupsP; }

  std::optional<Hdf5Location> commonLocation(const Hdf5Location& other) const {
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

  YAML::Node serialize() {
    YAML::Node node;
    node["file"] = fileP;
    node["group"] = groupsP;
    return node;
  }

  private:
  std::string fileP;
  std::vector<std::string> groupsP;
};

struct Hdf5AttributeWrite : public WriteInstruction {
  Hdf5Location location;
  std::string name;
  std::shared_ptr<writer::DataSource> dataSource;

  YAML::Node serialize() override {
    YAML::Node node;
    node["name"] = name;
    node["source"] = dataSource->serialize();
    node["location"] = location.serialize();
    node["writer"] = "hdf5";
    node["type"] = "attribute";
    return node;
  }

  Hdf5AttributeWrite(const Hdf5Location& location,
                     const std::string& name,
                     std::shared_ptr<writer::DataSource> dataSource)
      : location(location), name(name), dataSource(dataSource) {}

  explicit Hdf5AttributeWrite(YAML::Node node)
      : name(node["name"].as<std::string>()),
        dataSource(writer::DataSource::deserialize(node["source"])),
        location(Hdf5Location(node["location"])) {}

  std::vector<std::shared_ptr<DataSource>> dataSources() override { return {dataSource}; }
};

struct Hdf5DataWrite : public WriteInstruction {
  Hdf5Location location;
  std::string name;
  std::shared_ptr<writer::DataSource> dataSource;
  std::shared_ptr<datatype::Datatype> targetType;
  int compress;

  Hdf5DataWrite(const Hdf5Location& location,
                const std::string& name,
                std::shared_ptr<writer::DataSource> dataSource,
                std::shared_ptr<datatype::Datatype> targetType,
                int compress = 0)
      : location(location), name(name), dataSource(dataSource), targetType(targetType),
        compress(compress) {}

  YAML::Node serialize() override {
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

  explicit Hdf5DataWrite(YAML::Node node)
      : name(node["name"].as<std::string>()),
        dataSource(writer::DataSource::deserialize(node["source"])),
        location(Hdf5Location(node["location"])),
        targetType(datatype::Datatype::deserialize(node["targetType"])),
        compress(node["compress"].as<int>()) {}

  std::vector<std::shared_ptr<DataSource>> dataSources() override { return {dataSource}; }
};
} // namespace seissol::io::writer::instructions

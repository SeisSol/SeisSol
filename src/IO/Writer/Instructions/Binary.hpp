#pragma once

#include "Data.hpp"
#include "Instruction.hpp"
#include <memory>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {

struct BinaryWrite : public WriteInstruction {
  std::string filename;
  std::shared_ptr<writer::DataSource> dataSource;

  YAML::Node serialize() override {
    YAML::Node node;
    node["file"] = filename;
    node["source"] = dataSource->serialize();
    node["writer"] = "binary";
    return node;
  }

  BinaryWrite(const std::string& filename, std::shared_ptr<writer::DataSource> dataSource)
      : filename(filename), dataSource(dataSource) {}

  explicit BinaryWrite(YAML::Node node)
      : filename(node["file"].as<std::string>()),
        dataSource(writer::DataSource::deserialize(node["source"])) {}

  std::vector<std::shared_ptr<DataSource>> dataSources() override { return {dataSource}; }
};

} // namespace seissol::io::writer::instructions

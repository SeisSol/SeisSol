#include "Binary.hpp"

#include "Data.hpp"
#include "Instruction.hpp"
#include <memory>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {
YAML::Node BinaryWrite::serialize() {
  YAML::Node node;
  node["file"] = filename;
  node["source"] = dataSource->serialize();
  node["writer"] = "binary";
  return node;
}

BinaryWrite::BinaryWrite(const std::string& filename,
                         std::shared_ptr<writer::DataSource> dataSource)
    : filename(filename), dataSource(dataSource) {}

BinaryWrite::BinaryWrite(YAML::Node node)
    : filename(node["file"].as<std::string>()),
      dataSource(writer::DataSource::deserialize(node["source"])) {}

std::vector<std::shared_ptr<DataSource>> BinaryWrite::dataSources() { return {dataSource}; }

} // namespace seissol::io::writer::instructions

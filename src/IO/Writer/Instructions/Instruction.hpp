#pragma once

#include <IO/Writer/Instructions/Data.hpp>
#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {

class WriteInstruction {
  public:
  virtual YAML::Node serialize() = 0;
  static std::shared_ptr<WriteInstruction> deserialize(YAML::Node node);

  virtual std::vector<std::shared_ptr<DataSource>> dataSources() = 0;
};

} // namespace seissol::io::writer::instructions

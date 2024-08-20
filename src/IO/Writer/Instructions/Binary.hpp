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

  YAML::Node serialize() override;

  BinaryWrite(const std::string& filename, std::shared_ptr<writer::DataSource> dataSource);

  explicit BinaryWrite(YAML::Node node);

  std::vector<std::shared_ptr<DataSource>> dataSources() override;
};

} // namespace seissol::io::writer::instructions

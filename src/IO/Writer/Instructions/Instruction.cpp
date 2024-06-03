#include "Instruction.hpp"
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <memory>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {

std::shared_ptr<WriteInstruction> WriteInstruction::deserialize(YAML::Node node) {
  if (node["writer"].as<std::string>() == "hdf5") {
    if (node["type"].as<std::string>() == "data") {
      return std::make_shared<Hdf5DataWrite>(node);
    }
    if (node["type"].as<std::string>() == "attribute") {
      return std::make_shared<Hdf5AttributeWrite>(node);
    }
  }
  return std::shared_ptr<WriteInstruction>(nullptr);
}

} // namespace seissol::io::writer::instructions

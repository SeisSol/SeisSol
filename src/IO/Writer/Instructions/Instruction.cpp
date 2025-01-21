// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Instruction.h"
#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Hdf5.h>
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
  if (node["writer"].as<std::string>() == "binary") {
    return std::make_shared<BinaryWrite>(node);
  }
  return std::shared_ptr<WriteInstruction>(nullptr);
}

} // namespace seissol::io::writer::instructions

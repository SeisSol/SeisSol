// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_INSTRUCTION_HPP_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_INSTRUCTION_HPP_

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

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_INSTRUCTION_HPP_

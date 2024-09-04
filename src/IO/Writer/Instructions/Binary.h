// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_BINARY_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_BINARY_H_

#include "Data.h"
#include "Instruction.h"
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

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_BINARY_H_

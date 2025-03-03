// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_INSTRUCTION_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_INSTRUCTION_H_

#include <IO/Writer/Instructions/Data.h>
#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {

class WriteInstruction {
  public:
  virtual ~WriteInstruction() = default;
  virtual YAML::Node serialize() = 0;
  static std::shared_ptr<WriteInstruction> deserialize(YAML::Node node);

  virtual std::vector<std::shared_ptr<DataSource>> dataSources() = 0;
};

} // namespace seissol::io::writer::instructions

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_INSTRUCTION_H_

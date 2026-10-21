// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_BINARY_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_BINARY_H_

#include "Data.h"
#include "Instruction.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer::instructions {

struct BinaryWrite : public WriteInstruction {
  ~BinaryWrite() override = default;
  std::string filename;
  std::shared_ptr<writer::DataSource> dataSource;
  std::size_t alignment;
  bool append;
  //! Where in the file data that is not distributed goes; the file stays as it is around it.
  //! Without it, the data goes behind what the file holds.
  std::optional<std::size_t> offset;

  YAML::Node serialize() override;

  BinaryWrite(const std::string& filename,
              std::shared_ptr<writer::DataSource> dataSource,
              std::size_t alignment,
              bool append,
              std::optional<std::size_t> offset = std::nullopt);

  explicit BinaryWrite(YAML::Node node);

  std::vector<std::shared_ptr<DataSource>> dataSources() override;
};

} // namespace seissol::io::writer::instructions

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_BINARY_H_

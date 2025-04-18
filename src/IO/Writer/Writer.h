// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_WRITER_H_
#define SEISSOL_SRC_IO_WRITER_WRITER_H_

#include "Instructions/Instruction.h"
#include <IO/Writer/File/BinaryWriter.h>
#include <IO/Writer/File/Hdf5Writer.h>
#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <memory>
#include <yaml-cpp/yaml.h>

namespace async {
class ExecInfo;
} // namespace async

namespace seissol::io::writer {

class WriteInstance {
  public:
  WriteInstance(MPI_Comm comm);

  void write(const async::ExecInfo& info,
             const std::shared_ptr<instructions::WriteInstruction>& instruction);

  void close();

  private:
  file::Hdf5Writer hdf5;
  file::BinaryWriter binary;
};

class Writer {
  public:
  Writer();

  explicit Writer(const std::string& data);

  void addInstruction(const std::shared_ptr<instructions::WriteInstruction>& instruction);

  std::string serialize();

  WriteInstance beginWrite(const async::ExecInfo& info);

  void endWrite();

  [[nodiscard]] const std::vector<std::shared_ptr<instructions::WriteInstruction>>&
      getInstructions() const;

  private:
  std::vector<std::shared_ptr<instructions::WriteInstruction>> instructions;
};

struct ScheduledWriter {
  std::string name;
  double interval;
  std::function<Writer(const std::string&, std::size_t, double)> planWrite;
};

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_WRITER_H_

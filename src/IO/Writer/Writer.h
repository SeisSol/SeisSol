// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_WRITER_WRITER_HPP_
#define SEISSOL_SRC_IO_WRITER_WRITER_HPP_

#include "Instructions/Instruction.hpp"
#include "async/ExecInfo.h"
#include <IO/Writer/File/BinaryWriter.hpp>
#include <IO/Writer/File/Hdf5Writer.hpp>
#include <IO/Writer/Instructions/Binary.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

class WriteInstance {
  public:
  WriteInstance(MPI_Comm comm);

  void write(const async::ExecInfo& info,
             std::shared_ptr<instructions::WriteInstruction> instruction);

  void close();

  private:
  file::Hdf5Writer hdf5;
  file::BinaryWriter binary;
};

class Writer {
  public:
  Writer();

  explicit Writer(const std::string& data);

  void addInstruction(std::shared_ptr<instructions::WriteInstruction> instruction);

  std::string serialize();

  WriteInstance beginWrite(const async::ExecInfo& info);

  void endWrite();

  const std::vector<std::shared_ptr<instructions::WriteInstruction>>& getInstructions() const;

  private:
  std::vector<std::shared_ptr<instructions::WriteInstruction>> instructions;
};

struct ScheduledWriter {
  std::string name;
  double interval;
  std::function<Writer(const std::string&, std::size_t, double)> planWrite;
};

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_WRITER_HPP_

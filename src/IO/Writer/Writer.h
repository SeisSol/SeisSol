// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_WRITER_H_
#define SEISSOL_SRC_IO_WRITER_WRITER_H_

#include "IO/Writer/File/BinaryWriter.h"
#include "IO/Writer/File/Hdf5Writer.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "Instructions/Instruction.h"

#include <memory>
#include <mpi.h>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace async {
class ExecInfo;
} // namespace async

namespace seissol::io::writer {

/**
 * @brief Carries out the write plans of one output, one after the other.
 *
 * The files every plan writes stay open from one plan to the next (see file::Hdf5Writer::endWrite);
 * close() closes them.
 */
class WriteInstance {
  public:
  explicit WriteInstance(MPI_Comm comm, file::RunFiles* runFiles = nullptr);

  void write(const async::ExecInfo& info,
             const std::vector<std::shared_ptr<instructions::WriteInstruction>>& plan);

  void close();

  private:
  void write(const async::ExecInfo& info,
             const std::shared_ptr<instructions::WriteInstruction>& instruction);

  file::Hdf5Writer hdf5_;
  file::BinaryWriter binary_;
};

class Writer {
  public:
  Writer();

  explicit Writer(const std::string& data);

  void addInstruction(const std::shared_ptr<instructions::WriteInstruction>& instruction);
  void addInstructions(
      const std::vector<std::shared_ptr<instructions::WriteInstruction>>& instructions) {
    for (const auto& instruction : instructions) {
      addInstruction(instruction);
    }
  }

  std::string serialize();

  //! Carries out the plan with @p instance , which carried out the plans before it.
  void write(const async::ExecInfo& info, WriteInstance& instance) const;

  [[nodiscard]] const std::vector<std::shared_ptr<instructions::WriteInstruction>>&
      getInstructions() const;

  //! @brief The HDF5 files the plan writes to, each once, in the order they first appear in.
  [[nodiscard]] std::vector<std::string> hdf5Files() const;

  private:
  std::vector<std::shared_ptr<instructions::WriteInstruction>> instructions_;
};

struct ScheduledWriter {
  std::string name;
  double interval{};
  std::function<Writer(const std::string&, std::size_t, double)> planWrite;
};

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_WRITER_H_

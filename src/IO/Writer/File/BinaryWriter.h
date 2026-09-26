// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_
#define SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_
#include "IO/Writer/File/RunFiles.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"

#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::io::writer::file {

class BinaryFile {
  public:
  explicit BinaryFile(MPI_Comm comm);
  /**
   * @brief Opens @p name for writing; unless @p append , from its start.
   *
   * With @p backUp , an existing file is kept under a backup name first.
   */
  void openFile(const std::string& name, bool append, bool backUp = false);
  void writeGlobal(const void* data, std::size_t size);
  void writeDistributed(const void* data, std::size_t size);
  void align(std::size_t alignment);
  void closeFile();

  private:
  MPI_Comm comm_{MPI_COMM_NULL};
  MPI_File file_{MPI_FILE_NULL};
};

class BinaryWriter {
  public:
  /**
   * @brief A writer for one write plan; @p runFiles as for Hdf5Writer.
   */
  explicit BinaryWriter(MPI_Comm comm, RunFiles* runFiles = nullptr);

  void write(const async::ExecInfo& info, const instructions::BinaryWrite& write);

  void finalize();

  private:
  std::unordered_map<std::string, std::unique_ptr<BinaryFile>> openFiles_;
  MPI_Comm comm_{MPI_COMM_NULL};
  RunFiles* runFiles_{nullptr};
};
} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_

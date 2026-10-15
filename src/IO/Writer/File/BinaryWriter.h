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
#include <map>
#include <memory>
#include <mpi.h>
#include <set>
#include <stack>
#include <string>
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
  void writeGlobalAt(const void* data, std::size_t size, std::size_t offset);
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
   * @brief A writer for the write plans of one output; @p runFiles as for Hdf5Writer.
   */
  explicit BinaryWriter(MPI_Comm comm, RunFiles* runFiles = nullptr);

  void write(const async::ExecInfo& info, const instructions::BinaryWrite& write);

  //! Ends the writes of one plan, keeping the files open that every plan writes (cf. Hdf5Writer).
  void endWrite();

  //! Closes every file that is still open.
  void finalize();

  private:
  //! By name, so that every rank closes them in the same order: closing is collective.
  std::map<std::string, std::unique_ptr<BinaryFile>> openFiles_;
  //! The files the current plan wrote, and those the plan before wrote.
  std::set<std::string> written_;
  std::set<std::string> writtenBefore_;
  MPI_Comm comm_{MPI_COMM_NULL};
  RunFiles* runFiles_{nullptr};
};
} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_

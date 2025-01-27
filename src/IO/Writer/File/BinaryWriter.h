// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_
#define SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_
#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Data.h>
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
  BinaryFile(MPI_Comm comm);
  void openFile(const std::string& name);
  void writeGlobal(const void* data, std::size_t size);
  void writeDistributed(const void* data, std::size_t size);
  void closeFile();

  private:
  MPI_Comm comm;
  MPI_File file;
};

class BinaryWriter {
  public:
  BinaryWriter(MPI_Comm comm);

  void write(const async::ExecInfo& info, const instructions::BinaryWrite& write);

  void finalize();

  private:
  std::unordered_map<std::string, std::unique_ptr<BinaryFile>> openFiles;
  MPI_Comm comm;
};
} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_BINARYWRITER_H_

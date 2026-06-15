// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BinaryWriter.h"

#include "IO/Writer/Instructions/Binary.h"

#include <async/ExecInfo.h>
#include <cstddef>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <string>
#include <unordered_map>

namespace seissol::io::writer::file {

BinaryFile::BinaryFile(MPI_Comm comm) : comm_(comm) {}
void BinaryFile::openFile(const std::string& name) {
  MPI_File_open(comm_,
                name.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_APPEND,
                MPI_INFO_NULL,
                &file_);
}
void BinaryFile::writeGlobal(const void* data, std::size_t size) {
  int rank = 0;
  MPI_Comm_rank(comm_, &rank);
  if (rank == 0) {
    MPI_File_write(file_, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(comm_);
}
void BinaryFile::writeDistributed(const void* data, std::size_t size) {
  // TODO: get max write size
  MPI_File_write_all(file_, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
}
void BinaryFile::closeFile() { MPI_File_close(&file_); }

BinaryWriter::BinaryWriter(MPI_Comm comm) : comm_(comm) {}

void BinaryWriter::write(const async::ExecInfo& info, const instructions::BinaryWrite& write) {
  if (openFiles_.find(write.filename) == openFiles_.end()) {
    openFiles_[write.filename] = std::make_unique<BinaryFile>(BinaryFile(comm_));
    openFiles_[write.filename]->openFile(write.filename);
  }

  const void* dataPointer = write.dataSource->getPointer(info);

  // TODO: add dimensions
  const auto dataSize = write.dataSource->count(info) * write.dataSource->datatype()->size();

  if (write.dataSource->distributed()) {
    openFiles_[write.filename]->writeDistributed(dataPointer, dataSize);
  } else {
    openFiles_[write.filename]->writeGlobal(dataPointer, dataSize);
  }
}

void BinaryWriter::finalize() {
  for (auto& [_, file] : openFiles_) {
    file->closeFile();
  }
}

} // namespace seissol::io::writer::file

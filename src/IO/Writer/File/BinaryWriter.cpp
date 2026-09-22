// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "BinaryWriter.h"

#include "FileProperties.h"
#include "IO/Writer/Instructions/Binary.h"
#include "Parallel/MPI.h"
#include "RunFiles.h"

#include <async/ExecInfo.h>
#include <cstddef>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <string>
#include <unordered_map>

namespace seissol::io::writer::file {

BinaryFile::BinaryFile(MPI_Comm comm) : comm_(comm) {}
void BinaryFile::openFile(const std::string& name, bool append, bool backUp) {
  if (backUp) {
    int rank = 0;
    MPI_Comm_rank(comm_, &rank);
    if (rank == 0) {
      backUpFile(name);
    }
    // no rank may open the file before rank 0 moved the earlier one out of the way
    MPI_Barrier(comm_);
  }
  const auto mode = append ? MPI_MODE_APPEND : 0;
  // the same hints the HDF5 backend uses; the payload of an Xdmf output goes through here
  MPI_File_open(
      comm_, name.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY | mode, outputMpioHints(), &file_);
  if (!append) {
    // not appending means writing the file from its start, so whatever it held before goes;
    // distributed writes are placed behind the current end of the file
    MPI_File_set_size(file_, 0);
  }
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
  // TODO: handle size > usable
  MPI_Offset filesize = 0;
  MPI_File_get_size(file_, &filesize);
  std::size_t offset = 0;
  std::size_t total = 0;
  MPI_Exscan(&size, &offset, 1, Mpi::castToMpiType<std::size_t>(), MPI_SUM, comm_);
  MPI_Allreduce(&size, &total, 1, Mpi::castToMpiType<std::size_t>(), MPI_SUM, comm_);
  offset += filesize;

  MPI_File_write_at_all(file_, offset, data, size, MPI_BYTE, MPI_STATUS_IGNORE);
  MPI_File_seek(file_, 0, MPI_SEEK_END);
}
void BinaryFile::align(std::size_t alignment) {
  MPI_Offset position = 0;
  MPI_File_get_position(file_, &position);

  const auto alignedPosition = ((position + alignment - 1) / alignment) * alignment;
  MPI_File_seek(file_, alignedPosition, MPI_SEEK_SET);
}
void BinaryFile::closeFile() { MPI_File_close(&file_); }

BinaryWriter::BinaryWriter(MPI_Comm comm, RunFiles* runFiles) : comm_(comm), runFiles_(runFiles) {}

void BinaryWriter::write(const async::ExecInfo& info, const instructions::BinaryWrite& write) {
  if (openFiles_.find(write.filename) == openFiles_.end()) {
    openFiles_[write.filename] = std::make_unique<BinaryFile>(BinaryFile(comm_));
    // as for HDF5: the first write of a resuming run keeps what an earlier run left as a backup
    const bool backUp =
        runFiles_ != nullptr && runFiles_->firstWrite(write.filename) && runFiles_->resumed;
    openFiles_[write.filename]->openFile(write.filename, write.append, backUp);
  }

  const void* dataPointer = write.dataSource->getPointer(info);

  // TODO: add dimensions
  const auto dataSize = write.dataSource->count(info) * write.dataSource->datatype()->size();

  if (write.alignment > 0) {
    openFiles_[write.filename]->align(write.alignment);
  }

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

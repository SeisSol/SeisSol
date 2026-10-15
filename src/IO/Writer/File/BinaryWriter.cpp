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

#include <algorithm>
#include <async/ExecInfo.h>
#include <cstddef>
#include <hdf5.h>
#include <map>
#include <memory>
#include <mpi.h>
#include <set>
#include <string>
#include <utility>
#include <utils/logger.h>

namespace seissol::io::writer::file {

BinaryFile::BinaryFile(MPI_Comm comm, bool rootOnly, std::size_t maxPiece)
    : comm_(comm), rootOnly_(rootOnly), maxPiece_(std::max<std::size_t>(1, maxPiece)) {
  MPI_Comm_rank(comm_, &rank_);
}
bool BinaryFile::rootOnly() const { return rootOnly_; }
void BinaryFile::openFile(const std::string& name, bool append, bool backUp) {
  if (rootOnly_ && rank_ != 0) {
    return;
  }
  if (backUp) {
    if (rank_ == 0) {
      backUpFile(name);
    }
    // no rank may open the file before rank 0 moved the earlier one out of the way
    if (!rootOnly_) {
      MPI_Barrier(comm_);
    }
  }
  const auto mode = append ? MPI_MODE_APPEND : 0;
  // the same hints the HDF5 backend uses; the payload of an Xdmf output goes through here
  MPI_File_open(rootOnly_ ? MPI_COMM_SELF : comm_,
                name.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY | mode,
                outputMpioHints(),
                &file_);
  if (!append) {
    // not appending means writing the file from its start, so whatever it held before goes;
    // distributed writes are placed behind the current end of the file
    MPI_File_set_size(file_, 0);
  }
}
void BinaryFile::writeGlobal(const void* data, std::size_t size) {
  if (rank_ == 0) {
    const auto* bytes = static_cast<const char*>(data);
    for (std::size_t written = 0; written < size;) {
      const auto piece = std::min(size - written, maxPiece_);
      MPI_File_write(file_, bytes + written, static_cast<int>(piece), MPI_BYTE, MPI_STATUS_IGNORE);
      written += piece;
    }
  }
  // a distributed write behind it takes the size of the file for where it starts
  if (!rootOnly_) {
    MPI_Barrier(comm_);
  }
}
void BinaryFile::writeGlobalAt(const void* data, std::size_t size, std::size_t offset) {
  if (rank_ == 0) {
    const auto* bytes = static_cast<const char*>(data);
    for (std::size_t written = 0; written < size;) {
      const auto piece = std::min(size - written, maxPiece_);
      MPI_File_write_at(file_,
                        static_cast<MPI_Offset>(offset) + static_cast<MPI_Offset>(written),
                        bytes + written,
                        static_cast<int>(piece),
                        MPI_BYTE,
                        MPI_STATUS_IGNORE);
      written += piece;
    }
  }
  if (!rootOnly_) {
    MPI_Barrier(comm_);
  }
}
void BinaryFile::writeDistributed(const void* data, std::size_t size) {
  if (rootOnly_) {
    logError() << "Binary writer error: distributed data was written to a file that rank 0 writes "
                  "alone.";
  }
  MPI_Offset filesize = 0;
  MPI_File_get_size(file_, &filesize);
  std::size_t offset = 0;
  MPI_Exscan(&size, &offset, 1, Mpi::castToMpiType<std::size_t>(), MPI_SUM, comm_);
  offset += filesize;

  // A rank writes its data in pieces of at most maxPiece_ bytes (see BinaryFile::MaxPieceBytes).
  // The write is collective, so every rank takes part in as many rounds as the rank with the most
  // pieces, with nothing once its own are written.
  std::size_t rounds = (size + maxPiece_ - 1) / maxPiece_;
  MPI_Allreduce(MPI_IN_PLACE, &rounds, 1, Mpi::castToMpiType<std::size_t>(), MPI_MAX, comm_);
  const auto* bytes = static_cast<const char*>(data);
  for (std::size_t round = 0; round < rounds; ++round) {
    const auto start = std::min(size, round * maxPiece_);
    const auto piece = std::min(size - start, maxPiece_);
    MPI_File_write_at_all(file_,
                          static_cast<MPI_Offset>(offset) + static_cast<MPI_Offset>(start),
                          bytes + start,
                          static_cast<int>(piece),
                          MPI_BYTE,
                          MPI_STATUS_IGNORE);
  }
  MPI_File_seek(file_, 0, MPI_SEEK_END);
}
void BinaryFile::align(std::size_t alignment) {
  if (rootOnly_ && rank_ != 0) {
    return;
  }
  MPI_Offset position = 0;
  MPI_File_get_position(file_, &position);

  const auto alignedPosition = ((position + alignment - 1) / alignment) * alignment;
  MPI_File_seek(file_, alignedPosition, MPI_SEEK_SET);
}
void BinaryFile::closeFile() {
  if (rootOnly_ && rank_ != 0) {
    return;
  }
  MPI_File_close(&file_);
}

BinaryWriter::BinaryWriter(MPI_Comm comm, RunFiles* runFiles) : comm_(comm), runFiles_(runFiles) {}

void BinaryWriter::beginWrite(const std::set<std::string>& distributedFiles) {
  distributedFiles_ = distributedFiles;
}

void BinaryWriter::write(const async::ExecInfo& info, const instructions::BinaryWrite& write) {
  written_.insert(write.filename);
  if (openFiles_.find(write.filename) == openFiles_.end()) {
    const bool rootOnly = distributedFiles_.count(write.filename) == 0;
    openFiles_[write.filename] = std::make_unique<BinaryFile>(comm_, rootOnly);
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

  if (write.offset.has_value()) {
    if (write.dataSource->distributed()) {
      logError() << "Binary writer error: distributed data cannot be written to an offset in"
                 << write.filename;
    }
    openFiles_[write.filename]->writeGlobalAt(dataPointer, dataSize, write.offset.value());
  } else if (write.dataSource->distributed()) {
    openFiles_[write.filename]->writeDistributed(dataPointer, dataSize);
  } else {
    openFiles_[write.filename]->writeGlobal(dataPointer, dataSize);
  }
}

void BinaryWriter::endWrite() {
  for (auto entry = openFiles_.begin(); entry != openFiles_.end();) {
    const auto& name = entry->first;
    if (written_.count(name) > 0 && writtenBefore_.count(name) > 0) {
      ++entry;
    } else {
      entry->second->closeFile();
      entry = openFiles_.erase(entry);
    }
  }
  writtenBefore_ = std::move(written_);
  written_.clear();
}

void BinaryWriter::finalize() {
  for (auto& [_, file] : openFiles_) {
    file->closeFile();
  }
  openFiles_.clear();
  written_.clear();
  writtenBefore_.clear();
}

} // namespace seissol::io::writer::file

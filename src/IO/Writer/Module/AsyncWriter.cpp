// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AsyncWriter.h"

#include "IO/Writer/File/RunFiles.h"
#include "IO/Writer/Writer.h"

#include <async/ExecInfo.h>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <mutex>
#include <string>
#include <utils/logger.h>

namespace seissol::io::writer::module {
std::mutex AsyncWriter::globalLock = std::mutex();
std::condition_variable AsyncWriter::turn = std::condition_variable();
std::uint64_t AsyncWriter::nextTicket = 1;
file::RunFiles AsyncWriter::runFiles = file::RunFiles();

AsyncWriter::AsyncWriter() = default;

AsyncWriter::~AsyncWriter() = default;

void AsyncWriter::setName(const std::string& name) { name_ = name; }

void AsyncWriter::execInit(const async::ExecInfo& info, const AsyncWriterInit& params) {
  // (do nothing here)
}
void AsyncWriter::exec(const async::ExecInfo& info, const AsyncWriterExec& params) {
  {
    // for the Hdf5 implementations, we'll need to serialize writes
    // (TODO: make one AsyncWriter only in total)
    std::unique_lock lock(globalLock);
    if (params.ticket != 0) {
      turn.wait(lock, [&] { return nextTicket == params.ticket; });
    }
    if (params.close) {
      if (instance_.has_value()) {
        instance_.value().close();
        instance_.reset();
      }
    } else {
      const auto start = std::chrono::steady_clock::now();
      const void* data = info.buffer(PlanId);
      const size_t size = info.bufferSize(PlanId);
      const char* strData = reinterpret_cast<const char*>(data);

      // the ranks that write: all of them, or the executors of the nodes
      MPI_Comm comm = info.comm();

      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      if (printPlan_ && rank == 0) {
        logInfo() << "Printing current plan:" << std::string(strData, strData + size);
      }

      writer_ = Writer(std::string(strData, strData + size));
      runFiles.resumed = params.resumed;
      // one instance for all writes of the output, which keeps the files open that every write
      // goes to
      if (!instance_.has_value()) {
        instance_.emplace(comm, &runFiles);
      }
      writer_.write(info, instance_.value());
      // what the write itself took, on the thread or process that carries it out
      logInfo() << "Output Writer" << name_ << ": write done in"
                << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                << "s";
    }
    if (params.ticket != 0) {
      ++nextTicket;
    }
  }
  turn.notify_all();
}
void AsyncWriter::execWait(const async::ExecInfo& info) {
  // TODO: async finalize
}
void AsyncWriter::closeFiles() {
  const std::unique_lock lock(globalLock);
  if (instance_.has_value()) {
    instance_.value().close();
    instance_.reset();
  }
}
} // namespace seissol::io::writer::module

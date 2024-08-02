#include "AsyncWriter.hpp"

#include "async/ExecInfo.h"
#include <IO/Writer/Writer.hpp>
#include <cstddef>
#include <mpi.h>
#include <mutex>
#include <string>
#include <utils/logger.h>

namespace {} // namespace

namespace seissol::io::writer::module {
std::mutex AsyncWriter::globalLock = std::mutex();

void AsyncWriter::execInit(const async::ExecInfo& info, const AsyncWriterInit& params) {
  // (do nothing here)
}
void AsyncWriter::exec(const async::ExecInfo& info, const AsyncWriterExec& params) {
  const void* data = info.buffer(planId);
  const size_t size = info.bufferSize(planId);
  const char* strData = reinterpret_cast<const char*>(data);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (printPlan && rank == 0) {
    logInfo(rank) << "Printing current plan:" << std::string(strData, strData + size);
  }

  {
    // for the Hdf5 implementations, we'll need to serialize writes
    // (TODO: make one AsyncWriter only in total)
    std::lock_guard lock(globalLock);
    writer = Writer(std::string(strData, strData + size));
    instance = std::optional(writer.beginWrite(info));
    // for now write synchronously
    instance.value().close();
    instance.reset();
  }
}
void AsyncWriter::execWait(const async::ExecInfo& info) {
  // TODO: async finalize
}
void AsyncWriter::finalize() {}
} // namespace seissol::io::writer::module

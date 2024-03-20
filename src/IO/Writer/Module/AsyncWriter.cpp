#include "AsyncWriter.hpp"

#include "async/Module.h"
#include "async/ExecInfo.h"
#include <IO/Writer/Writer.hpp>
#include <string>

namespace seissol::io::writer::module {
void AsyncWriter::execInit(const async::ExecInfo& info, const AsyncWriterInit& params) {
  // (do nothing here)
}
void AsyncWriter::exec(const async::ExecInfo& info, const AsyncWriterExec& params) {
  const void* data = info.buffer(PlanId);
  size_t size = info.bufferSize(PlanId);
  const char* strData = reinterpret_cast<const char*>(data);

  //  AsyncIO::groupComm()
  logInfo() << info.buffer(0);
  logInfo() << info.buffer(1);
  logInfo() << info.buffer(2);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (printPlan && rank == 0) {
    logInfo(rank) << "Printing current plan:" << std::string(strData, strData + size);
  }
  writer = Writer(std::string(strData, strData + size));
  instance = std::optional(writer.beginWrite(info));
  // for now write synchronously
  instance.value().close();
  instance.reset();
}
void AsyncWriter::execWait(const async::ExecInfo& info, const AsyncWriterExec& params) {
  // TODO: async finalize
}
void AsyncWriter::finalize() {}
} // namespace seissol::io::writer::module

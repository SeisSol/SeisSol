// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "AsyncWriter.h"

#include "IO/Writer/Writer.h"

#include <async/ExecInfo.h>
#include <cstddef>
#include <mpi.h>
#include <mutex>
#include <string>
#include <utils/logger.h>

namespace seissol::io::writer::module {
std::mutex AsyncWriter::globalLock = std::mutex();

AsyncWriter::AsyncWriter() = default;

AsyncWriter::~AsyncWriter() = default;

void AsyncWriter::setComm(MPI_Comm comm) {
  // TODO: currently non-functional for ASYNC_MODE=mpi
  MPI_Comm_dup(comm, &this->comm);
}

void AsyncWriter::execInit(const async::ExecInfo& info, const AsyncWriterInit& params) {
  // (do nothing here)
}
void AsyncWriter::exec(const async::ExecInfo& info, const AsyncWriterExec& /*params*/) {
  const void* data = info.buffer(PlanId);
  const size_t size = info.bufferSize(PlanId);
  const char* strData = reinterpret_cast<const char*>(data);

  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (printPlan && rank == 0) {
    logInfo() << "Printing current plan:" << std::string(strData, strData + size);
  }

  {
    // for the Hdf5 implementations, we'll need to serialize writes
    // (TODO: make one AsyncWriter only in total)
    const std::scoped_lock lock(globalLock);
    writer = Writer(std::string(strData, strData + size));
    instance = std::optional(writer.beginWrite(info, comm));
    // for now write synchronously
    instance.value().close();
    instance.reset();
  }
}
void AsyncWriter::execWait(const async::ExecInfo& info) {
  // TODO: async finalize
}
void AsyncWriter::finalize() {
  if (comm != MPI_COMM_WORLD) {
    // comm has been duplicated
    MPI_Comm_free(&comm);
    comm = MPI_COMM_WORLD;
  }
}
} // namespace seissol::io::writer::module

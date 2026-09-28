// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_MODULE_ASYNCWRITER_H_
#define SEISSOL_SRC_IO_WRITER_MODULE_ASYNCWRITER_H_

#include "IO/Writer/File/Hdf5Writer.h"
#include "IO/Writer/Writer.h"

#include <async/ExecInfo.h>
#include <async/Module.h>
#include <condition_variable>
#include <cstdint>
#include <mutex>

namespace seissol::io::writer::module {
struct AsyncWriterInit {};

struct AsyncWriterExec {
  //! Whether the run resumes from a checkpoint (see file::RunFiles).
  bool resumed{false};
  /**
   * @brief The place of this write among the writes of all outputs, or zero for none.
   *
   * With one executor thread per output, the writes run in this order. Every write is collective
   * on the communicator of its output, and all of them share one lock; taking it in the order the
   * threads happen to arrive in could let two ranks each wait in a different write for the other.
   * The places are handed out on the main thread, in the same order on every rank.
   */
  std::uint64_t ticket{0};
};

class AsyncWriter {
  public:
  AsyncWriter();
  ~AsyncWriter();

  AsyncWriter(const AsyncWriter&) = delete;
  AsyncWriter(AsyncWriter&&) = delete;
  auto operator=(const AsyncWriter&) = delete;
  auto operator=(AsyncWriter&&) = delete;

  void setComm(MPI_Comm comm);
  void execInit(const async::ExecInfo& info, const AsyncWriterInit& params);
  void exec(const async::ExecInfo& info, const AsyncWriterExec& params);
  void execWait(const async::ExecInfo& info);
  void finalize();

  private:
  static constexpr int PlanId = 0;
  bool printPlan_{false};
  seissol::io::writer::Writer writer_;
  std::optional<seissol::io::writer::WriteInstance> instance_;
  MPI_Comm comm_{MPI_COMM_WORLD};

  static std::mutex globalLock;
  //! Signaled whenever a write is done, for the one whose turn is next.
  static std::condition_variable turn;
  //! The ticket of the write whose turn it is; guarded by globalLock.
  static std::uint64_t nextTicket;
  //! What the outputs of this run have written; shared by all of them, and guarded by globalLock.
  static file::RunFiles runFiles;
};

using AsyncWriterModule = async::Module<AsyncWriter, AsyncWriterInit, AsyncWriterExec>;
} // namespace seissol::io::writer::module

#endif // SEISSOL_SRC_IO_WRITER_MODULE_ASYNCWRITER_H_

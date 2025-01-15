// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_MODULE_ASYNCWRITER_H_
#define SEISSOL_SRC_IO_WRITER_MODULE_ASYNCWRITER_H_

#include "async/ExecInfo.h"
#include "async/Module.h"
#include <IO/Writer/Writer.h>
#include <mutex>

namespace seissol::io::writer::module {
struct AsyncWriterInit {};

struct AsyncWriterExec {};

class AsyncWriter {
  public:
  AsyncWriter() = default;
  void execInit(const async::ExecInfo& info, const AsyncWriterInit& params);
  void exec(const async::ExecInfo& info, const AsyncWriterExec& params);
  void execWait(const async::ExecInfo& info);
  void finalize();

  private:
  static constexpr int PlanId = 0;
  bool printPlan{false};
  seissol::io::writer::Writer writer;
  std::optional<seissol::io::writer::WriteInstance> instance;

  static std::mutex globalLock;
};

using AsyncWriterModule = async::Module<AsyncWriter, AsyncWriterInit, AsyncWriterExec>;
} // namespace seissol::io::writer::module

#endif // SEISSOL_SRC_IO_WRITER_MODULE_ASYNCWRITER_H_

#pragma once

#include "async/ExecInfo.h"
#include "async/Module.h"
#include <IO/Writer/Writer.hpp>
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
  static constexpr int planId = 0;
  bool printPlan{false};
  seissol::io::writer::Writer writer;
  std::optional<seissol::io::writer::WriteInstance> instance;

  static std::mutex globalLock;
};

using AsyncWriterModule = async::Module<AsyncWriter, AsyncWriterInit, AsyncWriterExec>;
} // namespace seissol::io::writer::module

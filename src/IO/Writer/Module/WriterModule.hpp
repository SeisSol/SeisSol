#pragma once

#include "AsyncWriter.hpp"
#include "Modules/Module.h"
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

namespace seissol::io::writer::module {

struct BufferPointer {
  size_t size;
  int id;
};

class WriterModule : public seissol::Module, private AsyncWriterModule {
  public:
  WriterModule(const ScheduledWriter& settings);
  void startup();
  void setUp() override;
  void simulationStart() override;
  void syncPoint(double time) override;
  void simulationEnd() override;
  void shutdown() override;

  private:
  unsigned planId;
  AsyncWriter executor;
  std::unordered_map<void*, BufferPointer> pointerMap;
  std::unordered_map<std::size_t, std::vector<int>> bufferMap;
  ScheduledWriter settings;
  double lastWrite;
};

} // namespace seissol::io::writer::module

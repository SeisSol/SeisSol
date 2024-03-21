#pragma once

#include "Writer/Module/WriterModule.hpp"
#include <IO/Instance/Checkpoint/CheckpointManager.hpp>
#include <IO/Writer/Writer.hpp>
#include <Modules/Module.h>
#include <Modules/Modules.h>
#include <memory>
#include <vector>

namespace seissol::io {

class OutputManager : public seissol::Module {
  public:
  OutputManager() = default;

  void setup() { Modules::registerHook(*this, ModuleHook::PostMPIInit); }

  void postMPIInit() override {}

  void addOutput(const writer::ScheduledWriter& writer) {
    modules.emplace_back(std::make_unique<seissol::io::writer::module::WriterModule>(writer));
    modules.back()->startup();
  }

  // TODO: de-couple checkpoint loading/writing (forward the CheckpointManager)
  
  void loadCheckpoint(const std::string& path) {
    checkpointManager.loadCheckpoint(path);
  }

  void setupCheckpoint(const std::string& path, double interval) {
    writer::ScheduledWriter checkpointScheduled;
    checkpointScheduled.name = "checkpoint";
    checkpointScheduled.planWrite = checkpointManager.makeWriter();
    checkpointScheduled.interval = interval;
    addOutput(checkpointScheduled);
  }

  private:
  instance::checkpoint::CheckpointManager checkpointManager;
  std::vector<std::unique_ptr<seissol::io::writer::module::WriterModule>> modules;
};

} // namespace seissol::io

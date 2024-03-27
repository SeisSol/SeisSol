#pragma once

#include "Writer/Module/WriterModule.hpp"
#include <IO/Instance/Checkpoint/CheckpointManager.hpp>
#include <IO/Writer/Writer.hpp>
#include <Modules/Module.h>
#include <Modules/Modules.h>
#include <memory>
#include <vector>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::io {

class OutputManager : public seissol::Module {
  public:
  OutputManager(SeisSol& seissolInstance);

  void setup() { Modules::registerHook(*this, ModuleHook::PostMPIInit); }

  void postMPIInit() override {}

  void addOutput(const writer::ScheduledWriter& writer);

  // TODO: de-couple checkpoint loading/writing (forward the CheckpointManager)

  void loadCheckpoint(const std::string& path) { checkpointManager.loadCheckpoint(path); }

  void setupCheckpoint(const std::string& path, double interval) {
    writer::ScheduledWriter checkpointScheduled;
    checkpointScheduled.name = "checkpoint";
    checkpointScheduled.planWrite = checkpointManager.makeWriter();
    checkpointScheduled.interval = interval;
    addOutput(checkpointScheduled);
  }

  instance::checkpoint::CheckpointManager& getCheckpointManager() { return checkpointManager; }

  private:
  instance::checkpoint::CheckpointManager checkpointManager;
  std::vector<std::unique_ptr<seissol::io::writer::module::WriterModule>> modules;
  SeisSol& seissolInstance;
};

} // namespace seissol::io

#include "Manager.hpp"

#include "Writer/Module/WriterModule.hpp"
#include <IO/Writer/Writer.hpp>
#include <SeisSol.h>
#include <memory>
#include <vector>

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::io {

OutputManager::OutputManager(SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

void OutputManager::addOutput(const writer::ScheduledWriter& writer) {
  modules.emplace_back(std::make_unique<seissol::io::writer::module::WriterModule>(
      seissolInstance.getSeisSolParameters().output.prefix, writer, seissolInstance.getPinning()));
  modules.back()->startup();
}

void OutputManager::postMPIInit() {
  // init ASYNC here
}

void OutputManager::shutdown() {
  // uninit ASYNC here
}

double OutputManager::loadCheckpoint(const std::string& path) {
  return checkpointManager.loadCheckpoint(path);
}

void OutputManager::setupCheckpoint(double interval) {
  writer::ScheduledWriter checkpointScheduled;
  checkpointScheduled.name = "checkpoint";
  checkpointScheduled.planWrite = checkpointManager.makeWriter();
  checkpointScheduled.interval = interval;
  addOutput(checkpointScheduled);
}

instance::checkpoint::CheckpointManager& OutputManager::getCheckpointManager() {
  return checkpointManager;
}

} // namespace seissol::io

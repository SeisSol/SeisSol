// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Manager.h"

#include "Writer/Module/WriterModule.h"
#include <IO/Writer/Writer.h>
#include <SeisSol.h>
#include <memory>
#include <string>
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

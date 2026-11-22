// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Manager.h"

#include "Core/IO/OutputContext.h"
#include "Core/IO/Writer/Writer.h"
#include "Writer/Module/WriterModule.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::io {

OutputManager::OutputManager(OutputContext context) : context_(context) {}

void OutputManager::addOutput(const writer::ScheduledWriter& writer) {
  modules_.emplace_back(std::make_unique<seissol::io::writer::module::WriterModule>(
      context_.outputPrefix, writer, context_));
  modules_.back()->startup();
}

double OutputManager::loadCheckpoint(const std::string& path) {
  return checkpointManager_.loadCheckpoint(path);
}

void OutputManager::setupCheckpoint(double interval) {
  writer::ScheduledWriter checkpointScheduled;
  checkpointScheduled.name = "checkpoint";
  checkpointScheduled.planWrite = checkpointManager_.makeWriter();
  checkpointScheduled.interval = interval;
  addOutput(checkpointScheduled);
}

instance::checkpoint::CheckpointManager& OutputManager::getCheckpointManager() {
  return checkpointManager_;
}

} // namespace seissol::io

// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_MANAGER_H_
#define SEISSOL_SRC_IO_MANAGER_H_

#include "Writer/Module/WriterModule.h"
#include <IO/Instance/Checkpoint/CheckpointManager.h>
#include <IO/Writer/Writer.h>
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

  void addOutput(const writer::ScheduledWriter& writer);

  // TODO: de-couple checkpoint loading/writing (forward the CheckpointManager)

  double loadCheckpoint(const std::string& path);

  void setupCheckpoint(double interval);

  instance::checkpoint::CheckpointManager& getCheckpointManager();

  private:
  instance::checkpoint::CheckpointManager checkpointManager;
  std::vector<std::unique_ptr<seissol::io::writer::module::WriterModule>> modules;
  SeisSol& seissolInstance;
};

} // namespace seissol::io

#endif // SEISSOL_SRC_IO_MANAGER_H_

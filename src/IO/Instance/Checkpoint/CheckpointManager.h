// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_CHECKPOINT_CHECKPOINTMANAGER_H_
#define SEISSOL_SRC_IO_INSTANCE_CHECKPOINT_CHECKPOINTMANAGER_H_

#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <IO/Reader/Distribution.h>
#include <IO/Reader/File/Hdf5Reader.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Writer.h>
#include <Initializer/Tree/LTSTree.h>
#include <Initializer/Tree/Layer.h>
#include <string>
#include <unordered_map>

#include "utils/logger.h"

namespace seissol::io::instance::checkpoint {

struct CheckpointVariable {
  std::string name;
  void* data;
  std::shared_ptr<datatype::Datatype> datatype;
};

struct CheckpointTree {
  std::string name;
  initializer::LTSTree* tree;
  std::vector<std::size_t> ids;
  std::vector<CheckpointVariable> variables;
};

class CheckpointManager {
  public:
  void registerTree(const std::string& name,
                    initializer::LTSTree* tree,
                    const std::vector<std::size_t>& ids) {
    dataRegistry[tree].name = name;
    dataRegistry[tree].tree = tree;
    dataRegistry[tree].ids = ids;
  }

  template <typename T>
  void registerData(const std::string& name,
                    initializer::LTSTree* tree,
                    initializer::Variable<T> var) {
    if (var.mask != initializer::LayerMask(Ghost)) {
      logError() << "Invalid layer mask for a checkpointing variable (i.e.: NYI).";
    }
    dataRegistry[tree].variables.emplace_back(
        CheckpointVariable{name, tree->var(var), datatype::inferDatatype<T>()});
  }

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter();

  double loadCheckpoint(const std::string& file);

  private:
  std::unordered_map<initializer::LTSTree*, CheckpointTree> dataRegistry;
};

} // namespace seissol::io::instance::checkpoint

#endif // SEISSOL_SRC_IO_INSTANCE_CHECKPOINT_CHECKPOINTMANAGER_H_

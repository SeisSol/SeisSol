#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/Inference.hpp>
#include <IO/Datatype/MPIType.hpp>
#include <IO/Reader/Distribution.hpp>
#include <IO/Reader/File/Hdf5Reader.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <IO/Writer/Writer.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Layer.hpp>
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

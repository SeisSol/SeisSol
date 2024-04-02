#pragma once

#include <Geometry/MeshReader.h>
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
  template <typename F>
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

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() {
    auto dataRegistry = this->dataRegistry;
    return [=](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
      writer::Writer writer;
      const auto filename = std::string("checkpoint-") + std::to_string(counter) + ".h5";
      for (auto& [_, ckpTree] : dataRegistry) {
        const std::size_t cells = ckpTree.tree->getNumberOfCells(Ghost);
        assert(cells == ckpTree.ids.size());
        std::size_t totalCells;
        MPI_Allreduce(&cells,
                      &totalCells,
                      1,
                      datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                      MPI_SUM,
                      MPI::mpi.comm());
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {"checkpoint", ckpTree.name}),
            "__count",
            writer::WriteInline::create(totalCells)));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, {"checkpoint", ckpTree.name}),
            "__ids",
            writer::WriteBuffer::create(ckpTree.ids.data(), ckpTree.ids.size()),
            datatype::inferDatatype<std::size_t>()));
        for (auto& variable : ckpTree.variables) {
          writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, {"checkpoint", ckpTree.name}),
              variable.name,
              std::make_shared<writer::WriteBuffer>(
                  variable.data, cells, variable.datatype, std::vector<std::size_t>()),
              variable.datatype));
        }
      }
      return writer;
    };
  }

  void loadCheckpoint(const std::string& file) {
    std::size_t storesize = 1;
    void* datastore = std::malloc(1);

    auto reader = reader::file::Hdf5Reader(seissol::MPI::mpi.comm());
    reader.openFile(file);
    for (auto& [_, ckpTree] : dataRegistry) {
      reader.openGroup(ckpTree.name);
      auto distributor = reader::Distributor(seissol::MPI::mpi.comm());
      std::size_t totalCount = reader.readAttributeScalar<std::size_t>("__count");

      auto groupIds = reader.readData<std::size_t>("__ids");
      distributor.setup(totalCount, groupIds, ckpTree.ids);
      for (auto& variable : ckpTree.variables) {
        const std::size_t count = reader.dataCount(variable.name);
        std::size_t currsize = count * variable.datatype->size();
        if (currsize > storesize) {
          datastore = std::realloc(datastore, currsize);
          storesize = currsize;
        }
        reader.readDataRaw(datastore, variable.name, count, variable.datatype);
        distributor.distribute(variable.data, datastore, datatype::convertToMPI(variable.datatype));
      }

      reader.closeGroup();
    }
    reader.closeFile();

    std::free(datastore);
  }

  private:
  std::unordered_map<initializer::LTSTree*, CheckpointTree> dataRegistry;
  geometry::MeshReader* meshReader;
};

} // namespace seissol::io::instance::checkpoint

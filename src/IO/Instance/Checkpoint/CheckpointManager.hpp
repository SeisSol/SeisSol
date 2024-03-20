#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/Inference.hpp>
#include <IO/Reader/Distribution.hpp>
#include <IO/Writer/Writer.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Layer.hpp>
#include <unordered_map>
namespace seissol::io::instance::checkpoint {

struct CheckpointVariable {
  std::string name;
  void* data;
  std::shared_ptr<datatype::Datatype> datatype;
};

struct CheckpointTree {
  std::string name;
  initializer::LTSTree* tree;
  const std::size_t* ids;
  std::vector<CheckpointVariable> variables;
};

class CheckpointManager {
  public:
  void registerTree(const std::string& name, initializer::LTSTree* tree, const std::size_t* ids) {
    dataRegistry[tree].name = name;
    dataRegistry[tree].tree = tree;
    dataRegistry[tree].ids = ids;
  }

  template <typename T>
  void registerData(const std::string& name,
                    initializer::LTSTree* tree,
                    initializer::Variable<T> var) {
    dataRegistry[tree].variables.emplace_back(
        CheckpointVariable{name, tree->var(var), datatype::inferDatatype<T>()});
  }

  std::function<writer::Writer(double)> makeWriter() {
    auto dataRegistry = this->dataRegistry;
    return [dataRegistry](double time) -> writer::Writer {
      writer::Writer writer;
      for (auto& [_, ckpTree] : dataRegistry) {
        auto cells = ckpTree->tree.getNumberOfCells(Ghost);
        MPI_Allreduce(&totalCells,
                      &cells,
                      1,
                      datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                      MPI_SUM,
                      comm);
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint", ckpTree.name}),
            "__count",
            writer::WriteInline::create(totalCells)));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint", ckpTree.name}),
            "__ids",
            writer::ElementTransformBuffer(),
            datatype::inferDatatype<std::size_t>() std::vector<std::size_t>()));
        for (auto& variable : ckpTree.variables) {
          // TODO(David): assert that we don't have ghost cells
          writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location("checkpoint.h5", {"checkpoint", ckpTree.name}),
              variable.name,
              writer::WriteBuffer(variable.data, cells, variable.datatype),
              variable.datatype,
              std::vector<std::size_t>()));
        }
      }
    };
  }

  void loadCheckpoint(const std::string& file) {
    auto reader = reader::Hdf5File(seissol::MPI::mpi.comm());
    reader.openFile(file);
    for (auto& [_, ckpTree] : dataRegistry) {
      reader.openGroup(ckpTree.name);
      auto distributor = reader::Distributor(seissol::MPI::mpi.comm());
      std::size_t totalCount =
          reader.readAttribute("__count", datatype::inferDatatype<std::size_t>());
      auto groupIds = reader.readData("__ids", datatype::inferDatatype<std::size_t>());
      distributor.setup(totalCount, groupIds.data(), ckpTree.ids);
      for (auto& variable : ckpTree.variables) {
        auto readData = reader.readData(variable.name, variable.datatype);
        distributor.distribute(variable.data, readData.data(), variable.datatype);
      }
      reader.closeGroup();
    }
    reader.closeFile();
  }

  private:
  std::unordered_map<initializer::LTSTree*, CheckpointTree> dataRegistry;
};

} // namespace seissol::io::instance::checkpoint

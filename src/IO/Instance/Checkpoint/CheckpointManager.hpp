#pragma once

#include <Geometry/MeshReader.h>
#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/Inference.hpp>
#include <IO/Reader/Distribution.hpp>
#include <IO/Reader/File/Hdf5Reader.hpp>
#include <IO/Writer/Writer.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Layer.hpp>
#include <string>
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
  template <typename F>
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
    auto meshReader = this->meshReader;
    return [=](double time) -> writer::Writer {
      writer::Writer writer;
      auto timestr = std::to_string(time);
      timestr.replace(timestr.find("."), 1, "-");
      const auto filename = std::string("checkpoint-") + timestr + ".h5";
      for (auto& [_, ckpTree] : dataRegistry) {
        const std::size_t cells = ckpTree.tree->getNumberOfCells(Ghost);
        assert(cells == meshReader->getElements().size());
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
            writer::WriteInline::create(totalCells, {})));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, {"checkpoint", ckpTree.name}),
            "__ids",
            writer::ElementBuffer<std::size_t>(cells,
                                               1,
                                               std::vector<std::size_t>(),
                                               [&](std::size_t* target, std::size_t index) {
                                                 target[0] =
                                                     meshReader->getElements()[index].globalId;
                                               }),
            datatype::inferDatatype<std::size_t>()));
        for (auto& variable : ckpTree.variables) {
          // TODO(David): assert that we don't have ghost cells
          writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, {"checkpoint", ckpTree.name}),
              variable.name,
              writer::WriteBuffer(variable.data, cells, variable.datatype, {}),
              variable.datatype,
              std::vector<std::size_t>()));
        }
      }
    };
  }

  void loadCheckpoint(const std::string& file) {
    auto reader = reader::file::Hdf5Reader(seissol::MPI::mpi.comm());
    reader.openFile(file);
    for (auto& [_, ckpTree] : dataRegistry) {
      reader.openGroup(ckpTree.name);
      auto distributor = reader::Distributor(seissol::MPI::mpi.comm());
      std::size_t totalCount = reader.readAttributeScalar<std::size_t>("__count");
      auto groupIds = reader.readData<std::size_t>("__ids", {});
      std::vector<std::size_t> currentGlobalIds(meshReader->getElements().size());
      for (std::size_t i = 0; i < currentGlobalIds.size(); ++i) {
        currentGlobalIds[i] = meshReader->getElements()[i].globalId;
      }
      distributor.setup(totalCount, groupIds.data(), currentGlobalIds.data());
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
  geometry::MeshReader* meshReader;
};

} // namespace seissol::io::instance::checkpoint

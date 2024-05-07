#include "CheckpointManager.hpp"

#include <IO/Datatype/Inference.hpp>
#include <IO/Datatype/MPIType.hpp>
#include <IO/Reader/Distribution.hpp>
#include <IO/Reader/File/Hdf5Reader.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <IO/Writer/Writer.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Layer.hpp>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <memory>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::io::instance::checkpoint {

std::function<writer::Writer(const std::string&, std::size_t, double)>
    CheckpointManager::makeWriter() {
  auto dataRegistry = this->dataRegistry;
  return [=](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
    writer::Writer writer;
    const auto filename = prefix + std::string("-checkpoint-") + std::to_string(counter) + ".h5";
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
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__time",
        writer::WriteInline::create(time)));
    return writer;
  };
}

double CheckpointManager::loadCheckpoint(const std::string& file) {
  std::size_t storesize = 1;
  void* datastore = std::malloc(1);

  auto reader = reader::file::Hdf5Reader(seissol::MPI::mpi.comm());
  reader.openFile(file);
  reader.openGroup("checkpoint");
  for (auto& [_, ckpTree] : dataRegistry) {
    reader.openGroup(ckpTree.name);
    auto distributor = reader::Distributor(seissol::MPI::mpi.comm());

    auto groupIds = reader.readData<std::size_t>("__ids");
    distributor.setup(groupIds, ckpTree.ids);
    for (auto& variable : ckpTree.variables) {
      const std::size_t count = reader.dataCount(variable.name);
      const std::size_t currsize = count * variable.datatype->size();
      if (currsize > storesize) {
        datastore = std::realloc(datastore, currsize);
        storesize = currsize;
      }
      reader.readDataRaw(datastore, variable.name, count, variable.datatype);
      distributor.distribute(variable.data, datastore, datatype::convertToMPI(variable.datatype));
    }

    reader.closeGroup();
  }
  const double time = reader.readAttributeScalar<double>("__time");
  reader.closeGroup();
  reader.closeFile();

  std::free(datastore);

  return time;
}

} // namespace seissol::io::instance::checkpoint

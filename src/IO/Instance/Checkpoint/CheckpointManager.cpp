// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CheckpointManager.h"

#include <Common/Constants.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <IO/Reader/Distribution.h>
#include <IO/Reader/File/Hdf5Reader.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Writer.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <Parallel/MPI.h>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "utils/logger.h"

namespace seissol::io::instance::checkpoint {

std::function<writer::Writer(const std::string&, std::size_t, double)>
    CheckpointManager::makeWriter() {
  auto dataRegistry = this->dataRegistry;
  return [=](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
    writer::Writer writer;
    const auto filename = prefix + std::string("-checkpoint-") + std::to_string(counter) + ".h5";
    for (const auto& [_, ckpTree] : dataRegistry) {
      const std::size_t cells = ckpTree.tree->getNumberOfCells(Ghost);
      assert(cells == ckpTree.ids.size());
      std::size_t totalCells = 0;
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
      for (const auto& variable : ckpTree.variables) {
        std::shared_ptr<writer::DataSource> dataSource;
        if (variable.pack.has_value()) {
          // data needs to be transformed, copy
          const auto elemSize = variable.memoryDatatype->size();
          const auto packFn = variable.pack.value();
          dataSource = writer::GeneratedBuffer::createElementwise<char>(
              elemSize,
              1,
              std::vector<std::size_t>(),
              [=](void* target, std::size_t index) {
                std::invoke(packFn,
                            target,
                            reinterpret_cast<const void*>(
                                reinterpret_cast<const char*>(variable.data) + index * elemSize));
              },
              variable.datatype);
        } else {
          // no transform; write-through
          dataSource = std::make_shared<writer::WriteBuffer>(
              variable.data, cells, variable.datatype, std::vector<std::size_t>());
        }
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, {"checkpoint", ckpTree.name}),
            variable.name,
            dataSource,
            variable.datatype));
      }
    }
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__time",
        writer::WriteInline::create(time)));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__order",
        writer::WriteInline::create(ConvergenceOrder)));
    return writer;
  };
}

double CheckpointManager::loadCheckpoint(const std::string& file) {
  std::size_t storesize = 0;
  void* datastore = nullptr;

  logInfo() << "Loading checkpoint...";
  logInfo() << "Checkpoint file:" << file;

  auto reader = reader::file::Hdf5Reader(seissol::MPI::mpi.comm());
  reader.openFile(file);
  reader.openGroup("checkpoint");
  const auto convergenceOrderRead = reader.readAttributeScalar<int>("__order");
  if (convergenceOrderRead != ConvergenceOrder) {
    logError() << "Convergence order does not match. Read:" << convergenceOrderRead;
  }
  for (auto& [_, ckpTree] : dataRegistry) {
    reader.openGroup(ckpTree.name);
    auto distributor = reader::Distributor(seissol::MPI::mpi.comm());

    logInfo() << "Reading group IDs for" << ckpTree.name;
    auto groupIds = reader.readData<std::size_t>("__ids");
    distributor.setup(groupIds, ckpTree.ids);

    std::vector<reader::Distributor::DistributionInstance> distributions;
    distributions.reserve(ckpTree.variables.size());
    for (auto& variable : ckpTree.variables) {
      logInfo() << "Reading variable" << ckpTree.name << "/" << variable.name;
      const std::size_t count = reader.dataCount(variable.name);
      const std::size_t currsize = count * variable.datatype->size();
      if (currsize > storesize) {
        std::free(datastore);
        datastore = std::malloc(currsize);
        if (datastore == nullptr) {
          logError() << "Realloc failed; maybe you are reading too much (checkpoint) data?";
        }
        storesize = currsize;

        // touch memory explicitly
        std::memset(datastore, 0, storesize);
      }
      reader.readDataRaw(datastore, variable.name, count, variable.datatype);
      const auto distribution =
          distributor.distributeRaw(variable.data,
                                    datastore,
                                    datatype::convertToMPI(variable.datatype),
                                    datatype::convertToMPI(variable.memoryDatatype),
                                    variable.unpack);
      distributions.push_back(distribution);
    }
    logInfo() << "Finishing data distribution for" << ckpTree.name;
    for (auto& distribution : distributions) {
      distribution.complete();
    }

    reader.closeGroup();
  }
  const auto time = reader.readAttributeScalar<double>("__time");
  reader.closeGroup();
  reader.closeFile();

  logInfo() << "Checkpoint loading complete.";

  std::free(datastore);

  return time;
}

} // namespace seissol::io::instance::checkpoint

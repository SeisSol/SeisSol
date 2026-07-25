// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CheckpointManager.h"

#include "Common/Constants.h"
#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Reader/Distribution.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "Parallel/MPI.h"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <mpi.h>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::instance::checkpoint {

std::function<writer::Writer(const std::string&, std::size_t, double)>
    CheckpointManager::makeWriter() {
  auto dataRegistry = this->dataRegistry_;
  return [=](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
    writer::Writer writer;
    const auto filename = prefix + std::string("-checkpoint-") + std::to_string(counter) + ".h5";

    for (const auto& [_, ckpStorage] : dataRegistry) {

      const auto location =
          writer::instructions::Hdf5Location(filename, {"checkpoint", ckpStorage.name});

      const std::size_t cells = ckpStorage.cells;
      assert(cells == ckpStorage.ids.size());
      std::size_t totalCells = 0;
      MPI_Allreduce(&cells,
                    &totalCells,
                    1,
                    datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                    MPI_SUM,
                    Mpi::mpi.comm());
      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          location,
          "__ids",
          writer::WriteBuffer::create(ckpStorage.ids.data(), ckpStorage.ids.size()),
          datatype::inferDatatype<std::size_t>()));
      for (const auto& variable : ckpStorage.variables) {
        std::shared_ptr<writer::DataSource> dataSource;
        if (variable.pack.has_value()) {
          // data needs to be transformed, copy
          const auto elemSize = variable.memoryDatatype->size();
          const auto packFn = variable.pack.value();
          dataSource = writer::GeneratedBuffer::createElementwise<std::byte>(
              elemSize,
              1,
              std::vector<std::size_t>(),
              [=](void* target, std::size_t index) {
                std::invoke(
                    packFn,
                    target,
                    reinterpret_cast<const void*>(
                        reinterpret_cast<const std::byte*>(variable.data) + index * elemSize));
              },
              variable.datatype);
        } else {
          // no transform; write-through
          dataSource = std::make_shared<writer::WriteBuffer>(
              variable.data, cells, variable.datatype, std::vector<std::size_t>());
        }
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            location, variable.name, dataSource, variable.datatype));
      }
    }

    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__order",
        writer::WriteInline::create(ConvergenceOrder)));

    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__time",
        writer::WriteInline::create(time)));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__alignment",
        writer::WriteInline::create(Alignment)));
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
        writer::instructions::Hdf5Location(filename, {"checkpoint"}),
        "__version",
        writer::WriteInline::create(FormatVersion)));

    return writer;
  };
}

double CheckpointManager::loadCheckpoint(const std::string& file) {
  std::size_t storesize = 0;
  void* datastore = nullptr;

  logInfo() << "Loading checkpoint...";
  logInfo() << "Checkpoint file:" << file;

  auto reader = reader::file::Hdf5Reader(seissol::Mpi::mpi.comm());
  reader.openFile(file);
  reader.openGroup("checkpoint");
  const auto versionRead = reader.readAttributeScalar<std::size_t>("__version");
  if (versionRead != FormatVersion) {
    logWarning() << "The checkpoint format version does not match. Read:" << versionRead
                 << "vs build:" << FormatVersion << ". The checkpoint read will most likely fail.";
  }
  const auto convergenceOrderRead = reader.readAttributeScalar<std::size_t>("__order");
  if (convergenceOrderRead != ConvergenceOrder) {
    logWarning() << "Convergence order does not match. Read:" << convergenceOrderRead
                 << "vs build:" << ConvergenceOrder
                 << ". The checkpoint read will most likely fail.";
  }
  const auto alignmentRead = reader.readAttributeScalar<std::size_t>("__alignment");
  if (alignmentRead != Alignment) {
    logWarning() << "Memory alignment and padding does not match. Read:" << alignmentRead
                 << "vs build:" << Alignment << ". The checkpoint read will most likely fail.";
  }

  // the distributors should not be moved (for the closures in the DistributionInstance not to point
  // to invalid memory).
  std::list<reader::Distributor> distributors;
  std::vector<reader::Distributor::DistributionInstance> distributions;

  for (auto& [_, ckpStorage] : dataRegistry_) {
    if (reader.hasEntry(ckpStorage.name)) {
      reader.openGroup(ckpStorage.name);

      logInfo() << "Reading group IDs for" << ckpStorage.name;
      auto groupIds = reader.readData<std::size_t>("__ids");
      auto& distributor = distributors.emplace_back(seissol::Mpi::mpi.comm());
      distributor.setup(groupIds, ckpStorage.ids);

      for (auto& variable : ckpStorage.variables) {
        if (reader.hasEntry(variable.name)) {
          logInfo() << "Reading variable" << ckpStorage.name << "/" << variable.name;
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
        } else {
          logWarning() << "The following variable was not found checkpoint file (instead, using "
                          "default values):"
                       << ckpStorage.name << "/" << variable.name;
        }
      }
      reader.closeGroup();
    } else {
      logWarning() << "The following data storage was not found in the checkpoint file:"
                   << ckpStorage.name;
    }
  }
  const auto time = reader.readAttributeScalar<double>("__time");
  reader.closeGroup();
  reader.closeFile();

  logInfo() << "Finishing data distribution.";
  for (auto& distribution : distributions) {
    distribution.complete();
  }

  logInfo() << "Checkpoint loading complete.";

  std::free(datastore);

  return time;
}

} // namespace seissol::io::instance::checkpoint

// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WriterModule.h"

#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Module/AsyncWriter.h"
#include "IO/Writer/Writer.h"
#include "Modules/Modules.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/Pin.h"
#include "SeisSol.h"

#include <cassert>
#include <cmath>
#include <cstring>
#include <optional>
#include <string>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::writer::module {

WriterModule::WriterModule(const std::string& prefix,
                           const ScheduledWriter& settings,
                           const parallel::Pinning& pinning,
                           SeisSol& seissolInstance)
    : rank_(seissol::Mpi::mpi.rank()), prefix_(prefix), settings_(settings), pinning_(pinning),
      seissolInstance_(seissolInstance) {}

void WriterModule::setUp() {
  logInfo() << "Output Writer" << settings_.name << ": setup.";
  executor_.setComm(seissol::Mpi::mpi.comm());
  setExecutor(executor_);
  // TODO: adjust the CommThread call here
  if (isAffinityNecessary() && useCommThread(seissol::Mpi::mpi, seissolInstance_.env())) {
    const auto freeCpus = pinning_.getFreeCPUsMask();
    logInfo() << "Output Writer" << settings_.name
              << ": thread affinity: " << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
    setAffinityIfNecessary(freeCpus);
  }
}

void WriterModule::startup() {
  logInfo() << "Output Writer" << settings_.name << ": startup, running at interval"
            << settings_.interval;
  init();

  // we want ASYNC to like us, hence we need to enter a non-zero size here
  planId_ = addBuffer(nullptr, 1, true);
  assert(planId_ == 0);

  callInit(AsyncWriterInit{});
  // TODO: pinning

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  Modules::registerHook(*this, ModuleHook::SimulationEnd);
  Modules::registerHook(*this, ModuleHook::Shutdown);
  setSyncInterval(settings_.interval);
}

void WriterModule::simulationStart(std::optional<double> checkpointTime) {
  if (checkpointTime.value_or(0) == 0) {
    syncPoint(0);
  }
}

void WriterModule::syncPoint(double time) {
  if (lastWrite_ >= 0) {
    logInfo() << "Output Writer" << settings_.name << ": finishing previous write from"
              << lastWrite_;
  }
  wait();
  logInfo() << "Output Writer" << settings_.name << ": preparing write at" << time;

  // request the write plan
  auto writeCount = static_cast<int>(std::round(time / syncInterval()));
  auto writer = settings_.planWrite(prefix_, writeCount, time);

  // prepare the data in the plan
  std::unordered_set<DataSource*> handledSources;
  std::vector<int> idsToSend;
  std::unordered_set<int> idSet;
  for (const auto& instruction : writer.getInstructions()) {
    for (auto& dataSource : instruction->dataSources()) {
      if (handledSources.find(dataSource.get()) == handledSources.end()) {
        // TODO: make a better flag than distributed here
        if (dataSource->distributed()) {
          // NOTE: the following structure is suboptimal, because it's not respecting basic OOP
          // practices. Not sure, if it's important to really care about that... But there may be
          // more beautiful ways for some day.
          auto id = [&]() -> int {
            if (dynamic_cast<WriteBuffer*>(dataSource.get()) != nullptr) {
              // pass-through buffer
              auto* writeBuffer = dynamic_cast<WriteBuffer*>(dataSource.get());
              const auto* pointer = writeBuffer->getLocalPointer();
              auto size = writeBuffer->getLocalSize();
              if (pointerMap_.find(pointer) == pointerMap_.end()) {
                BufferPointer repr;
                repr.id = addBuffer(pointer, size);
                repr.size = size;
                pointerMap_[pointer] = repr;
                return repr.id;
              } else {
                auto& repr = pointerMap_.at(pointer);
                if (repr.size != size) {
                  if (idSet.find(repr.id) != idSet.end()) {
                    // it is ok to request the same buffer multiple times, but not with different
                    // sizes (that's currently still unsupported)
                    logError() << "The same buffer is requested with different sizes. This is not "
                                  "supported at the moment.";
                  }
                  resizeBuffer(repr.id, pointer, size);
                  repr.size = size;
                }
                return repr.id;
              }
            }
            if (dynamic_cast<AdhocBuffer*>(dataSource.get()) != nullptr) {
              // managed buffer
              // for now, recycle existing buffers of the same size
              // this does of course assume that we don't change the size too often...
              auto* adhocBuffer = dynamic_cast<AdhocBuffer*>(dataSource.get());
              auto targetSize = adhocBuffer->getTargetSize();
              const auto foundId = [&]() -> int {
                for (auto id : bufferMap_[targetSize]) {
                  if (idSet.find(id) == idSet.end()) {
                    return id;
                  }
                }
                auto newId = addBuffer(nullptr, targetSize);
                bufferMap_[targetSize].push_back(newId);
                return newId;
              }();

              // avoid assert in ASYNC by checking targetSize == 0 explicitly
              void* bufferPtr = targetSize == 0 ? nullptr : managedBuffer<void*>(foundId);
              adhocBuffer->setData(bufferPtr);
              return foundId;
            }
            logError() << "Unsupported buffer type.";
            return -1;
          }();
          if (id == -1) {
            logError() << "Internal buffer error.";
          }
          idSet.emplace(id);
          idsToSend.push_back(id);
          dataSource->assignId(id);
          handledSources.emplace(dataSource.get());
        }
      }
    }
  }

  // take care of the plan (i.e. resize our managed buffer and send it)
  auto serialized = writer.serialize();
  resizeBuffer(planId_, nullptr, serialized.size());
  char* planBuffer = managedBuffer<char*>(planId_);
  std::memcpy(planBuffer, serialized.c_str(), serialized.size());
  sendBuffer(planId_, serialized.size());

  // send the plan data
  for (const int id : idsToSend) {
    sendBuffer(id);
  }

  logInfo() << "Output Writer" << settings_.name << ": triggering write at" << time;
  lastWrite_ = time;
  call(AsyncWriterExec{});
}

void WriterModule::simulationEnd() {
  logInfo() << "Output Writer" << settings_.name << ": finishing output";
  wait();
}

void WriterModule::shutdown() {
  logInfo() << "Output Writer" << settings_.name << ": shutdown";
  finalize();
}

} // namespace seissol::io::writer::module

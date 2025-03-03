// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WriterModule.h"
#include "utils/logger.h"
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Module/AsyncWriter.h>
#include <IO/Writer/Writer.h>
#include <Modules/Modules.h>
#include <Parallel/Helper.h>
#include <Parallel/MPI.h>
#include <Parallel/Pin.h>
#include <cassert>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_set>
#include <vector>

namespace seissol::io::writer::module {

WriterModule::WriterModule(const std::string& prefix,
                           const ScheduledWriter& settings,
                           const parallel::Pinning& pinning)
    : rank(seissol::MPI::mpi.rank()), prefix(prefix), settings(settings), pinning(pinning) {}

void WriterModule::setUp() {
  logInfo() << "Output Writer" << settings.name << ": setup.";
  setExecutor(executor);
  if (isAffinityNecessary() && useCommThread(seissol::MPI::mpi)) {
    const auto freeCpus = pinning.getFreeCPUsMask();
    logInfo() << "Output Writer" << settings.name
              << ": thread affinity: " << parallel::Pinning::maskToString(freeCpus);
    if (parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
    setAffinityIfNecessary(freeCpus);
  }
}

void WriterModule::startup() {
  logInfo() << "Output Writer" << settings.name << ": startup, running at interval"
            << settings.interval;
  init();

  // we want ASYNC to like us, hence we need to enter a non-zero size here
  planId = addBuffer(nullptr, 1, true);
  assert(planId == 0);

  callInit(AsyncWriterInit{});
  // TODO: pinning

  Modules::registerHook(*this, ModuleHook::SimulationStart);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
  Modules::registerHook(*this, ModuleHook::SimulationEnd);
  Modules::registerHook(*this, ModuleHook::Shutdown);
  setSyncInterval(settings.interval);
}

void WriterModule::simulationStart() { syncPoint(0); }

void WriterModule::syncPoint(double time) {
  if (lastWrite >= 0) {
    logInfo() << "Output Writer" << settings.name << ": finishing previous write from" << lastWrite;
  }
  wait();
  logInfo() << "Output Writer" << settings.name << ": preparing write at" << time;

  // request the write plan
  auto writeCount = static_cast<int>(std::round(time / syncInterval()));
  auto writer = settings.planWrite(prefix, writeCount, time);

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
              if (pointerMap.find(pointer) == pointerMap.end()) {
                BufferPointer repr;
                repr.id = addBuffer(pointer, size);
                repr.size = size;
                pointerMap[pointer] = repr;
                return repr.id;
              } else {
                auto& repr = pointerMap.at(pointer);
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
                for (auto id : bufferMap[targetSize]) {
                  if (idSet.find(id) == idSet.end()) {
                    return id;
                  }
                }
                auto newId = addBuffer(nullptr, targetSize);
                bufferMap[targetSize].push_back(newId);
                return newId;
              }();
              void* bufferPtr = managedBuffer<void*>(foundId);
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
  resizeBuffer(planId, nullptr, serialized.size());
  char* planBuffer = managedBuffer<char*>(planId);
  std::memcpy(planBuffer, serialized.c_str(), serialized.size());
  sendBuffer(planId, serialized.size());

  // send the plan data
  for (const int id : idsToSend) {
    sendBuffer(id);
  }

  logInfo() << "Output Writer" << settings.name << ": triggering write at" << time;
  lastWrite = time;
  call(AsyncWriterExec{});
}

void WriterModule::simulationEnd() {
  logInfo() << "Output Writer" << settings.name << ": finishing output";
  wait();
}

void WriterModule::shutdown() {
  logInfo() << "Output Writer" << settings.name << ": shutdown";
  finalize();
}

} // namespace seissol::io::writer::module

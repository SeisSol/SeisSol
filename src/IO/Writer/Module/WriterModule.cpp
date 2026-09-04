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
#include <cstring>
#include <optional>
#include <string>
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
              << ": thread affinity: " << parallel::Pinning::maskToString(freeCpus) << "("
              << parallel::Pinning::maskToStringShort(freeCpus).c_str() << ")";
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
  planId_ = AsyncWriterModule::addBuffer(nullptr, 1, true);
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
  auto writer = settings_.planWrite(prefix_, writeCount_, time);

  // decide which buffer each data source of the plan uses
  const auto idsToSend = registry_.assign(writer);

  // take care of the plan (i.e. resize our managed buffer and send it)
  auto serialized = writer.serialize();
  AsyncWriterModule::resizeBuffer(planId_, nullptr, serialized.size());
  auto* planBuffer = AsyncWriterModule::managedBuffer<char*>(planId_);
  std::memcpy(planBuffer, serialized.c_str(), serialized.size());
  sendBuffer(planId_, serialized.size());

  // send the plan data
  for (const int id : idsToSend) {
    sendBuffer(id);
  }

  logInfo() << "Output Writer" << settings_.name << ": triggering write at" << time;
  lastWrite_ = time;
  ++writeCount_;
  call(AsyncWriterExec{});
}

int WriterModule::addBuffer(const void* pointer, std::size_t size) {
  return AsyncWriterModule::addBuffer(pointer, size, pointer == nullptr);
}

void WriterModule::resizeBuffer(int id, const void* pointer, std::size_t size) {
  AsyncWriterModule::resizeBuffer(id, pointer, size);
}

void* WriterModule::managedBuffer(int id) { return AsyncWriterModule::managedBuffer<void*>(id); }

void WriterModule::simulationEnd() {
  logInfo() << "Output Writer" << settings_.name << ": finishing output";
  wait();
}

void WriterModule::shutdown() {
  logInfo() << "Output Writer" << settings_.name << ": shutdown";
  finalize();
}

} // namespace seissol::io::writer::module

// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "WriterModule.h"

#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Module/AsyncWriter.h"
#include "IO/Writer/Writer.h"
#include "Modules/Modules.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/Pin.h"
#include "SeisSol.h"

#include <async/Config.h>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mpi.h>
#include <optional>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::writer::module {

namespace {
/**
 * @brief The place of the next call among the writes of all outputs (see AsyncWriterExec::ticket).
 *
 * With one executor thread per output, on every rank or on the executor ranks of a node, the
 * writes of all outputs have to run in the order they are issued in on the main thread, which is
 * the same on every rank. Synchronously, they run in that order anyway.
 */
std::uint64_t nextTicket(async::Mode mode) {
  static std::uint64_t issued = 0;
  return mode == async::Mode::Sync ? 0 : ++issued;
}
} // namespace

WriterModule::WriterModule(std::string prefix,
                           ScheduledWriter settings,
                           const parallel::Pinning& pinning,
                           SeisSol& seissolInstance)
    : AsyncWriterModule(seissolInstance.asyncIO().config()), rank_(seissol::Mpi::mpi.rank()),
      prefix_(std::move(prefix)), settings_(std::move(settings)), pinning_(pinning),
      seissolInstance_(seissolInstance) {}

void WriterModule::setUp() {
  logInfo() << "Output Writer" << settings_.name << ": setup.";
  executor_.setName(settings_.name);
  setExecutor(executor_);
  // The thread goes to the CPUs no OpenMP worker uses, with a communication thread or without one:
  // left to ASYNC, it takes the last CPU of the mask it inherits from the main thread, which is the
  // main thread's own core, or its own CPU, while that core computes. The ranks of a node find
  // those CPUs together, so every rank looks them up, also where no thread of the output runs.
  const auto freeCpus = pinning_.getFreeCPUsMask();
  if (isAffinityNecessary()) {
    if (!parallel::Pinning::freeCPUsMaskEmpty(freeCpus)) {
      logInfo() << "Output Writer" << settings_.name
                << ": thread affinity: " << parallel::Pinning::maskToString(freeCpus) << "("
                << parallel::Pinning::maskToStringShort(freeCpus).c_str() << ")";
      setAffinityIfNecessary(freeCpus);
    } else if (useCommThread(seissol::Mpi::mpi, seissolInstance_.env())) {
      logError() << "There are no free CPUs left. Make sure to leave one for the I/O thread(s).";
    }
  }
}

void WriterModule::startup() {
  logInfo() << "Output Writer" << settings_.name << ": startup, running at interval"
            << settings_.interval;
  // the dispatcher runs already; the module joins it, together with all ranks
  seissolInstance_.asyncIO().attach(*this);
  init();

  // we want ASYNC to like us, hence we need to enter a non-zero size here
  planId_ = AsyncWriterModule::addBuffer(nullptr, 1, true);
  assert(planId_ == 0);

  callInit(AsyncWriterInit{});
  // TODO: pinning

  seissolInstance_.modules().registerHook(*this, ModuleHook::SimulationStart);
  seissolInstance_.modules().registerHook(*this, ModuleHook::SynchronizationPoint);
  seissolInstance_.modules().registerHook(*this, ModuleHook::SimulationEnd);
  seissolInstance_.modules().registerHook(*this, ModuleHook::Shutdown);
  setSyncInterval(settings_.interval);
}

void WriterModule::simulationStart(std::optional<double> checkpointTime) {
  const auto startTime = checkpointTime.value_or(0);
  resumed_ = startTime > 0;
  if (startTime == 0) {
    syncPoint(0);
  } else {
    // resume the numbering instead of starting over, which would overwrite the output of the run
    // this one continues
    writeCount_ = outputCountBefore(startTime, settings_.interval);
    logInfo() << "Output Writer" << settings_.name << ": resuming at output" << writeCount_;
  }
}

void WriterModule::syncPoint(double time) {
  const auto start = std::chrono::steady_clock::now();
  if (lastWrite_ >= 0) {
    logInfo() << "Output Writer" << settings_.name << ": finishing previous write from"
              << lastWrite_;
  }
  wait();
  const auto waited = std::chrono::steady_clock::now();
  logInfo() << "Output Writer" << settings_.name << ": preparing write at" << time;

  // request the write plan
  auto writer = settings_.planWrite(prefix_, writeCount_, time);

  // every HDF5 file of the output says which run wrote to it last
  if (!seissolInstance_.runId().empty()) {
    for (const auto& file : writer.hdf5Files()) {
      writer.addInstruction(std::make_shared<instructions::Hdf5AttributeWrite>(
          instructions::Hdf5Location(file, {}),
          "runId",
          WriteInline::createString(seissolInstance_.runId())));
    }
  }
  const auto planned = std::chrono::steady_clock::now();

  // decide which buffer each data source of the plan uses; a generated one fills its buffer here
  const auto idsToSend = registry_.assign(writer);
  const auto generated = std::chrono::steady_clock::now();

  // take care of the plan (i.e. resize our managed buffer and send it)
  const auto serialized = writer.serialize();
  AsyncWriterModule::resizeBuffer(planId_, nullptr, serialized.size());
  auto* planBuffer = AsyncWriterModule::managedBuffer<char*>(planId_);
  std::memcpy(planBuffer, serialized.c_str(), serialized.size());
  sendBuffer(planId_, serialized.size());

  // send the plan data
  for (const int id : idsToSend) {
    sendBuffer(id);
  }

  const auto handedOver = std::chrono::steady_clock::now();
  const auto seconds = [](auto from, auto to) {
    return std::chrono::duration<double>(to - from).count();
  };
  // where the main thread spent its time on this write, up to here
  logInfo() << "Output Writer" << settings_.name << ": triggering write at" << time << "after"
            << seconds(start, handedOver) << "s (waiting for the previous write"
            << seconds(start, waited) << "s, plan" << seconds(waited, planned) << "s, data"
            << seconds(planned, generated) << "s, handover" << seconds(generated, handedOver)
            << "s)";
  lastWrite_ = time;
  ++writeCount_;
  call(AsyncWriterExec{.resumed = resumed_,
                       .ticket = nextTicket(seissolInstance_.asyncIO().config().mode())});
}

int WriterModule::addBuffer(const void* pointer, std::size_t size) {
  // (every rank contributes its own data, also to a buffer ASYNC manages)
  return AsyncWriterModule::addBuffer(pointer, size);
}

void WriterModule::resizeBuffer(int id, const void* pointer, std::size_t size) {
  AsyncWriterModule::resizeBuffer(id, pointer, size);
}

void* WriterModule::managedBuffer(int id) { return AsyncWriterModule::managedBuffer<void*>(id); }

bool WriterModule::collective() const {
  // the executors of a node join the parts of its ranks under the same id
  return seissolInstance_.asyncIO().config().mode() == async::Mode::MPI;
}

void WriterModule::agreeOnAny(std::vector<int>& flags) {
  MPI_Allreduce(MPI_IN_PLACE,
                flags.data(),
                static_cast<int>(flags.size()),
                MPI_INT,
                MPI_LOR,
                seissol::Mpi::mpi.comm());
}

void WriterModule::simulationEnd() {
  logInfo() << "Output Writer" << settings_.name << ": finishing output";
  wait();
  // the files the writes kept open for the next one; closing is collective, and takes its turn
  // like a write
  call(AsyncWriterExec{.resumed = resumed_,
                       .ticket = nextTicket(seissolInstance_.asyncIO().config().mode()),
                       .close = true});
  wait();
}

void WriterModule::shutdown() {
  logInfo() << "Output Writer" << settings_.name << ": shutdown";
  finalize();
  // in case the simulation did not end the output
  executor_.closeFiles();
}

} // namespace seissol::io::writer::module

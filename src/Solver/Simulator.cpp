// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Sebastian Rettenberger

#include "Simulator.h"

#include "Kernels/Common.h"
#include "Memory/Tree/Layer.h"
#include "Modules/Modules.h"
#include "Modules/OutputWorker.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/RunSummary.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "Parallel/OpenMP.h"
#include "Parallel/Pin.h"
#include "ResultWriter/AnalysisWriter.h"
#include "ResultWriter/EnergyOutput.h"
#include "SeisSol.h"
#include "TimeStepping/TimeManager.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <mpi.h>
#include <optional>
#include <utils/logger.h>

namespace seissol {

Simulator::Simulator() = default;

void Simulator::setFinalTime(double finalTime) {
  assert(finalTime > 0);
  this->finalTime_ = finalTime;
}

void Simulator::setUsePlasticity(bool plasticity) { usePlasticity_ = plasticity; }

void Simulator::setCurrentTime(double currentTime) {
  assert(currentTime >= 0);
  this->currentTime_ = currentTime;
  checkpoint_ = true;
}

void Simulator::abort() { abortRequested_ = true; }

namespace {
/**
 * Lets the outputs process their synchronization points on a thread of their own, if asked for
 * and possible: the host has to keep copies of the data (see SEISSOL_USM), which the device does
 * not touch in the next epoch.
 */
void startOutputThread(SeisSol& seissolInstance) {
  if (!useOutputThread(seissolInstance.env())) {
    return;
  }
  if (!isDeviceOn() || useUSM() || !seissolInstance.timeManager().clustersOnDevice()) {
    logWarning() << "SEISSOL_OUTPUT_THREAD needs a device, all clusters on it, and separate buffers"
                 << "for the host (SEISSOL_USM=0); the outputs are processed on the main thread.";
    return;
  }
  // the CPUs the OpenMP threads of the main thread do not use, or else theirs, which are idle while
  // the device computes
  auto cpus = seissolInstance.pinning().getFreeCPUsMask();
  if (parallel::Pinning::freeCPUsMaskEmpty(cpus)) {
    cpus = parallel::Pinning::getWorkerUnionMask();
  }
  const auto threads =
      std::min(OutputWorker::cpuCount(cpus), static_cast<int>(OpenMP::threadCount()));
  Modules::startOutputWorker(cpus, threads);
  logInfo() << "The outputs are processed on a thread of their own, with" << threads
            << "OpenMP threads, on the CPUs" << parallel::Pinning::maskToString(cpus).c_str();
}
} // namespace

void Simulator::simulate(SeisSol& seissolInstance) {
  SCOREP_USER_REGION("simulate", SCOREP_USER_REGION_TYPE_FUNCTION)

  auto* faultOutputManager = seissolInstance.timeManager().faultOutputManager();
  faultOutputManager->writePickpointOutput(0.0, 0.0);

  auto& summary = seissolInstance.runSummary();
  summary["simulation"]["startTime"] = currentTime_;
  summary["simulation"]["endTime"] = finalTime_;

  Stopwatch simulationStopwatch;
  if (checkpoint_) {
    logInfo() << "The simulation will run from" << currentTime_ << "s (checkpoint time) until"
              << finalTime_ << "s.";
  } else {
    logInfo() << "The simulation will run until" << finalTime_ << "s.";
  }
  simulationStopwatch.start();

  Stopwatch computeStopwatch;
  Stopwatch ioStopwatch;

  ioStopwatch.start();

  // Set start time (required when loading a check)
  seissolInstance.timeManager().setInitialTimes(currentTime_);

  const double timeTolerance = seissolInstance.timeManager().getTimeTolerance();

  startOutputThread(seissolInstance);

  // Write initial wave field snapshot
  if (checkpoint_) {
    Modules::callSimulationStartHook(currentTime_);
  } else {
    Modules::callSimulationStartHook(std::optional<double>{});
  }

  // derive next synchronization time
  double upcomingTime = finalTime_;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min(upcomingTime, Modules::callSyncHook(currentTime_, 0.0));

  // The times of an epoch are the differences of what the stopwatches accumulated (pause returns
  // that); split measures from the last start, through any pause since, so after the IO of an
  // epoch the split of the compute stopwatch would hold the IO as well.
  double lastIo = ioStopwatch.pause();
  double lastCompute = 0;

  summary.setRankValue("simulationInitialIo", lastIo);
  summary["time"]["simulation"]["initialIo"] =
      monitoring::toJson(Stopwatch::print("Time spent for initial IO:", lastIo));

  double lastSplit = simulationStopwatch.split();

  // use an empty log message as visual separator
  logInfo() << "";

  while (finalTime_ > currentTime_ + timeTolerance) {
    if (upcomingTime < currentTime_ + timeTolerance) {
      logError() << "Simulator did not advance in time from" << currentTime_ << "s to"
                 << upcomingTime << "s.";
    }
    if (aborted_) {
      logInfo() << "Aborting simulation.";
      break;
    }

    // update the DOFs
    logInfo() << "Start simulation epoch. (from" << currentTime_ << "s to" << upcomingTime << "s)";
    const double epochStart = currentTime_;

    computeStopwatch.start();
    seissolInstance.timeManager().advanceInTime(upcomingTime);
    const double compute = computeStopwatch.pause();

    const auto rateDone = upcomingTime / finalTime_;
    const auto percentageDone = rateDone * 100;
    logInfo() << "End simulation epoch. (at" << upcomingTime << "s; total" << percentageDone
              << "% done)";

    ioStopwatch.start();

    // update current time
    currentTime_ = upcomingTime;

    // the outputs of the last synchronization point may still read the host copies of the data;
    // once they are done, what they decided (see abort()) is the same on all ranks
    Modules::waitForOutputs();
    aborted_ = aborted_ || abortRequested_;

    // the data the modules due now read on the host
    seissolInstance.timeManager().synchronizeToHost(
        Modules::hostDataAt(currentTime_, timeTolerance));

    // Check all synchronization point hooks and set the new upcoming time
    upcomingTime = std::min(finalTime_, Modules::callSyncHook(currentTime_, timeTolerance));
    if (!Modules::outputWorkerRunning()) {
      aborted_ = aborted_ || abortRequested_;
    }

    const double io = ioStopwatch.pause();

    const double currentSplit = simulationStopwatch.split();
    monitoring::EpochRecord epoch;
    epoch.timeFrom = epochStart;
    epoch.timeTo = currentTime_;
    epoch.wall = Stopwatch::print("Time spent this epoch (total):", currentSplit - lastSplit);
    epoch.compute = Stopwatch::print("Time spent this epoch (compute):", compute - lastCompute);
    epoch.io = Stopwatch::print("Time spent this epoch (blocking IO):", io - lastIo);
    lastCompute = compute;
    lastIo = io;

    const auto rates = seissolInstance.flopCounter().printPerformanceUpdate(currentSplit);
    epoch.hardwareFlopPerSecond = rates.hardwareFlopPerSecond;
    epoch.nonzeroFlopPerSecond = rates.nonzeroFlopPerSecond;
    epoch.kernelBytesPerSecond = rates.kernelBytesPerSecond;
    lastSplit = currentSplit;

    summary["time"]["simulation"]["elapsed"] =
        monitoring::toJson(Stopwatch::print("Total time spent:", currentSplit));
    Stopwatch::print("Expected time until finish:", currentSplit / rateDone * (1.0 - rateDone));
    summary.addEpoch(epoch, finalTime_);

    // use an empty log message as visual separator
    logInfo() << "";
  }

  // the outputs whose interval does not end at the final time write here
  ioStopwatch.start();

  // all data, for the outputs at the end and the analysis
  Modules::waitForOutputs();
  seissolInstance.timeManager().synchronizeTo(initializer::AllocationPlace::Host);

  Modules::callSyncHook(currentTime_, timeTolerance, true);
  // the outputs end on the main thread
  Modules::stopOutputWorker();

  // and the writes still under way finish
  Modules::callHook<ModuleHook::SimulationEnd>();

  const double ioTime = ioStopwatch.pause();

  const double wallTime = simulationStopwatch.pause();
  auto& times = summary["time"]["simulation"];
  summary.setRankValue("simulationTotal", wallTime);
  summary.setRankValue("simulationCompute", lastCompute);
  summary.setRankValue("simulationIo", ioTime);
  const auto totalTime = simulationStopwatch.printTime("Simulation time (total):");
  times["total"] = monitoring::toJson(totalTime);
  times["compute"] = monitoring::toJson(computeStopwatch.printTime("Simulation time (compute):"));
  times["io"] = monitoring::toJson(ioStopwatch.printTime("Simulation time (blocking IO):"));
  summary["simulation"]["currentTime"] = currentTime_;
  summary["simulation"]["aborted"] = aborted_;

  summary["kernels"] = seissolInstance.timeManager().printComputationTime();
  const double kernelTime = seissolInstance.timeManager().kernelTime();
  summary.setRankValue("kernels", kernelTime);
  // what the compute time holds besides the kernels: waiting for the other clusters and ranks, the
  // communication on the main thread, and the scheduling of the clusters
  times["outsideKernels"] = summary.overRanks("outsideKernels", lastCompute - kernelTime);

  seissolInstance.analysisWriter().printAnalysis(currentTime_);

  summary["performance"] = seissolInstance.flopCounter().printPerformanceSummary(wallTime);

  // the work the simulation did, over all ranks: cell updates, which compare across orders,
  // equations and machines where FLOP do not
  auto updates = seissolInstance.timeManager().updates();
  double localCellUpdates = 0;
  double localDuplicateUpdates = 0;
  double localFaceUpdates = 0;
  for (std::size_t cluster = 0; cluster < updates.cells.size(); ++cluster) {
    localCellUpdates += updates.cells[cluster];
    localDuplicateUpdates += updates.duplicateCells[cluster];
    localFaceUpdates += updates.dynamicRuptureFaces[cluster];
  }
  summary.setRankValue("cellUpdates", localCellUpdates);
  summary.setRankValue("duplicateCellUpdates", localDuplicateUpdates);
  summary.setRankValue("dynamicRuptureFaceUpdates", localFaceUpdates);
  for (auto* counts : {&updates.cells, &updates.duplicateCells, &updates.dynamicRuptureFaces}) {
    MPI_Allreduce(MPI_IN_PLACE,
                  counts->data(),
                  static_cast<int>(counts->size()),
                  MPI_DOUBLE,
                  MPI_SUM,
                  Mpi::mpi.comm());
  }
  double cellUpdates = 0;
  double duplicateUpdates = 0;
  double faceUpdates = 0;
  auto& clusters = summary["lts"]["clusters"];
  for (std::size_t cluster = 0; cluster < updates.cells.size(); ++cluster) {
    cellUpdates += updates.cells[cluster];
    duplicateUpdates += updates.duplicateCells[cluster];
    faceUpdates += updates.dynamicRuptureFaces[cluster];
    if (cluster < clusters.size()) {
      clusters.at(cluster)["cellUpdates"] = updates.cells[cluster];
      clusters.at(cluster)["duplicateCellUpdates"] = updates.duplicateCells[cluster];
      clusters.at(cluster)["dynamicRuptureFaceUpdates"] = updates.dynamicRuptureFaces[cluster];
    }
  }
  summary["simulation"]["cellUpdates"] = cellUpdates;
  summary["simulation"]["duplicateCellUpdates"] = duplicateUpdates;
  summary["simulation"]["dynamicRuptureFaceUpdates"] = faceUpdates;
  // the slowest rank decides how long the simulation takes
  summary["performance"]["cellUpdatesPerSecond"] = cellUpdates / totalTime.max;
  logInfo() << "Cell updates:" << cellUpdates << "(" << cellUpdates / totalTime.max
            << "per second); duplicate updates of copy cells:" << duplicateUpdates;
}

} // namespace seissol

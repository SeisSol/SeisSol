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
#include "Monitoring/DeviceTimeline.h"
#include "Monitoring/FlopCounter.h"
#include "Monitoring/Instrumentation.h"
#include "Monitoring/RunSummary.h"
#include "Monitoring/Stopwatch.h"
#include "Monitoring/Trace.h"
#include "Monitoring/Watchdog.h"
#include "Numerical/Statistics.h"
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
  // with unified memory, only the snapshots keep what the host reads apart from the device
  if (!isDeviceOn() || !seissolInstance.timeManager().clustersOnDevice() ||
      (useUSM() && !keepOutputSnapshots(seissolInstance.env()))) {
    logWarning()
        << "SEISSOL_OUTPUT_THREAD needs a device and all clusters on it, with SEISSOL_USM=1"
        << "also the snapshots (SEISSOL_OUTPUT_SNAPSHOT=1); the outputs are processed on"
        << "the main thread.";
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
  seissolInstance.modules().startOutputWorker(cpus, threads);
  logInfo() << "The outputs are processed on a thread of their own, with" << threads
            << "OpenMP threads, on the CPUs" << parallel::Pinning::maskToString(cpus).c_str();
}

/**
 * Prints how busy the threads and the device were in an epoch, and records it in the epoch.
 * Collective.
 */
void printActivity(const solver::EpochActivity& activity, monitoring::EpochRecord& epoch) {
  const double capacity = activity.threads * activity.wall;
  epoch.threadLoops =
      statistics::parallelSummary(capacity > 0 ? activity.threadLoops / capacity : 0.0);
  // with a device, the threads may run no loops at all
  if (!DeviceTimeline::active() || epoch.threadLoops.max > 0) {
    logInfo() << "Time the threads spent in loops this epoch:" << 100 * epoch.threadLoops.mean
              << "% (mean over all ranks; minimum" << 100 * epoch.threadLoops.min << "%)";
  }
  if (!DeviceTimeline::active()) {
    return;
  }
  const double span = activity.device.span;
  const auto share = [&](double time) {
    return statistics::parallelSummary(span > 0 ? time / span : 0.0);
  };
  epoch.deviceTimed = true;
  epoch.deviceComputing = share(activity.device.computing);
  epoch.deviceExchanges = share(activity.device.busyWith(DeviceTimeline::Kind::Exchanges));
  epoch.deviceIdle = share(span - activity.device.busy);
  logInfo() << "Time the device computed this epoch:" << 100 * epoch.deviceComputing.mean
            << "% (mean over all ranks; minimum" << 100 * epoch.deviceComputing.min
            << "%), exchanged the halo" << 100 * epoch.deviceExchanges.mean << "%, and was idle"
            << 100 * epoch.deviceIdle.mean << "%";
}
} // namespace

void Simulator::simulate(SeisSol& seissolInstance) {
  SCOREP_USER_REGION("simulate", SCOREP_USER_REGION_TYPE_FUNCTION)

  auto* faultOutputManager = seissolInstance.timeManager().faultOutputManager();
  faultOutputManager->recordInitialPickpointOutput(currentTime_);

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

  auto& modules = seissolInstance.modules();
  startOutputThread(seissolInstance);
  const bool snapshots =
      modules.outputWorkerRunning() && keepOutputSnapshots(seissolInstance.env());

  // Write initial wave field snapshot
  if (checkpoint_) {
    modules.callSimulationStartHook(currentTime_);
  } else {
    modules.callSimulationStartHook(std::optional<double>{});
  }

  // derive next synchronization time
  double upcomingTime = finalTime_;
  // NOTE: This will not call the module specific implementation of the synchronization hook
  // since the current time is the simulation start time. We only use this function here to
  // get correct upcoming time. To be on the safe side, we use zero time tolerance.
  upcomingTime = std::min(upcomingTime, modules.callSyncHook(currentTime_, 0.0));

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

  // the time stepping is where a trace starts (SEISSOL_TRACE)
  monitoring::Trace::start(monitoring::Trace::settingsFromEnv(seissolInstance.env()),
                           seissol::Mpi::mpi.comm());

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
    seissolInstance.timeManager().setWatchdogPhase(monitoring::Watchdog::Phase::Other,
                                                   "the outputs of the last epoch");
    modules.waitForOutputs();
    aborted_ = aborted_ || abortRequested_;
    seissolInstance.timeManager().setWatchdogPhase(monitoring::Watchdog::Phase::Other,
                                                   "a synchronization point");

    // the data the modules due now read on the host; where their thread processes the
    // synchronization point, it copies the data there from the snapshots, while the device goes on
    // with the next epoch
    const bool fromSnapshot =
        snapshots && modules.processesAsynchronouslyAt(currentTime_, timeTolerance);
    seissolInstance.timeManager().synchronizeToHost(modules.hostDataAt(currentTime_, timeTolerance),
                                                    fromSnapshot);

    // Check all synchronization point hooks and set the new upcoming time
    upcomingTime = std::min(
        finalTime_, modules.callSyncHook(currentTime_, timeTolerance, false, [&seissolInstance] {
          seissolInstance.timeManager().copySnapshotToHost();
        }));
    if (!modules.outputWorkerRunning()) {
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
    printActivity(seissolInstance.timeManager().lastEpochActivity(), epoch);

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
  seissolInstance.timeManager().setWatchdogPhase(monitoring::Watchdog::Phase::Other,
                                                 "the outputs at the end");
  modules.waitForOutputs();
  seissolInstance.timeManager().synchronizeTo(initializer::AllocationPlace::Host);

  modules.callSyncHook(currentTime_, timeTolerance, true);
  // the outputs end on the main thread
  modules.stopOutputWorker();

  // and the writes still under way finish
  modules.callHook<ModuleHook::SimulationEnd>();
  seissolInstance.timeManager().setWatchdogPhase(monitoring::Watchdog::Phase::Idle, "the end");

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

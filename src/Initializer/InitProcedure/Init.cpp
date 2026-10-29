// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Init.h"

#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Json.h"
#include "Geometry/MeshReader.h"
#include "InitIO.h"
#include "InitMesh.h"
#include "InitModel.h"
#include "InitSideConditions.h"
#include "Initializer/InitProcedure/InitLayout.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Memory/Tree/Layer.h"
#include "Monitoring/DeviceTimeline.h"
#include "Monitoring/Memory.h"
#include "Monitoring/RunSummary.h"
#include "Monitoring/Stopwatch.h"
#include "Monitoring/Trace.h"
#include "Monitoring/Unit.h"
#include "Numerical/Statistics.h"
#include "Parallel/Helper.h"
#include "Parallel/MPI.h"
#include "ResultWriter/ThreadsPinningWriter.h"
#include "SeisSol.h"

#include <async/Config.h>
#include <cstddef>
#include <string>
#include <utils/logger.h>

#ifdef ACL_DEVICE
#include <Device/device.h>
#include <ostream>
#include <sstream>
#endif

namespace seissol::initializer::initprocedure {

namespace {

void reportDeviceMemoryStatus() {
#ifdef ACL_DEVICE
  const device::DeviceInstance& device = device::DeviceInstance::instance();
  if (device.api().getCurrentlyOccupiedMem() > device.api().getMaxAvailableMem()) {
    std::stringstream stream;

    stream << "Memory of the device is overloaded." << std::endl
           << "Totally allocated device memory: "
           << UnitByte.formatPrefix(device.api().getCurrentlyOccupiedMem()) << std::endl
           << "Allocated unified memory: "
           << UnitByte.formatPrefix(device.api().getCurrentlyOccupiedUnifiedMem()) << std::endl
           << "Memory capacity of device: "
           << UnitByte.formatPrefix(device.api().getMaxAvailableMem());

    logError() << stream.str();
  } else {
    const double fraction = device.api().getCurrentlyOccupiedMem() /
                            static_cast<double>(device.api().getMaxAvailableMem());
    const auto summary = seissol::statistics::parallelSummary(fraction * 100.0);
    logInfo() << "occupied memory on devices (%):"
              << " mean =" << summary.mean << " std =" << summary.std << " min =" << summary.min
              << " median =" << summary.median << " max =" << summary.max;
  }
#endif
}

//! @brief How the cells are spread over the ranks. Collective.
json::Value meshSection(seissol::SeisSol& seissolInstance) {
  const auto& mesh = seissolInstance.meshReader();
  auto& summary = seissolInstance.runSummary();
  const auto rank = seissol::Mpi::mpi.rank();
  std::size_t copyCells = 0;
  for (const auto& element : mesh.getElements()) {
    if (geometry::isCopy(element, rank)) {
      ++copyCells;
    }
  }
  json::Value section;
  section["cells"] = summary.overRanks("cells", static_cast<double>(mesh.getElements().size()));
  section["copyCells"] = summary.overRanks("copyCells", static_cast<double>(copyCells));
  section["ghostCells"] =
      summary.overRanks("ghostCells", static_cast<double>(mesh.linearGhostlayer().size()));
  section["faultFaces"] =
      summary.overRanks("faultFaces", static_cast<double>(mesh.getFault().size()));
  section["neighborRanks"] =
      summary.overRanks("neighborRanks", static_cast<double>(mesh.getMPINeighbors().size()));
  return section;
}

void initSeisSol(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();
  auto& summary = seissolInstance.runSummary();

  // set g
  seissolInstance.gravitationSetup().acceleration = seissolParams.model.gravitationalAcceleration;

  // the stages that do not time themselves are timed here; the summary gets all of them
  const auto timed = [&](const char* name, const char* column, const char* message, auto&& stage) {
    seissol::Stopwatch watch;
    watch.start();
    stage(seissolInstance);
    summary.setRankValue(column, watch.pause());
    summary["time"]["init"][name] = monitoring::toJson(watch.printTime(message));
  };

  // initialization procedure
  summary.setStage("init.mesh");
  seissol::initializer::initprocedure::initMesh(seissolInstance);
  summary["mesh"] = meshSection(seissolInstance);
  summary.setStage("init.layout");
  timed("layout",
        "initLayout",
        "Layout initialized in:",
        seissol::initializer::initprocedure::initLayout);
  summary.setStage("init.model");
  seissol::initializer::initprocedure::initModel(seissolInstance);
  summary.setStage("init.sideConditions");
  timed("sideConditions",
        "initSideConditions",
        "Side conditions initialized in:",
        seissol::initializer::initprocedure::initSideConditions);
  summary.setStage("init.io");
  timed("io", "initIo", "Output initialized in:", seissol::initializer::initprocedure::initIO);

  // synchronize data to device
  seissolInstance.memoryManager().synchronizeTo(seissol::initializer::AllocationPlace::Device);

  // set up simulator
  auto& sim = seissolInstance.simulator();
  sim.setUsePlasticity(seissolParams.model.plasticity);
  sim.setFinalTime(seissolParams.timeStepping.endTime);
}

void reportHardwareRelatedStatus(seissol::SeisSol& seissolInstance) {
  reportDeviceMemoryStatus();

  const auto& seissolParams = seissolInstance.parameters();
  writer::ThreadsPinningWriter pinningWriter(seissolParams.output.prefix, seissolInstance.runId());
  pinningWriter.write(seissolInstance.pinning(), seissolInstance.env());

  auto& parallel = seissolInstance.runSummary()["parallel"];
  parallel["commThread"] = useCommThread(seissol::Mpi::mpi, seissolInstance.env());
  const auto layout = seissolInstance.asyncIO().layout();
  parallel["asyncMode"] = [&layout] {
    switch (layout.mode) {
    case async::Mode::Thread:
      return "thread";
    case async::Mode::MPI:
      return "node";
    default:
      return "sync";
    }
  }();
  parallel["asyncAuto"] = seissolInstance.asyncIO().automatic();
  if (layout.mode == async::Mode::MPI) {
    auto& executors = parallel["asyncExecutors"];
    executors["executors"] = layout.executors;
    executors["nodes"] = layout.nodes;
    executors["executorsPerNode"] = layout.executorsPerNode;
    executors["groupNodes"] = layout.groupNodes;
    executors["exclusive"] = layout.exclusive;
    executors["minGroupRanks"] = layout.minGroupRanks;
    executors["maxGroupRanks"] = layout.maxGroupRanks;
  }
}

/**
 * @brief The memory the run holds, per rank: the peak of the host memory, and what is allocated on
 * the device. Collective; prints the host peak as well.
 *
 * The table of the ranks gets the host peak and the occupied device memory, in the columns named
 * after @p moment.
 */
json::Value memorySection(monitoring::RunSummary& runSummary, const std::string& moment) {
  json::Value section;
  const auto peak = monitoring::peakHostMemory();
  if (peak.has_value()) {
    const auto localPeak = static_cast<double>(peak.value());
    runSummary.setRankValue("hostMemory" + moment, localPeak);
    const auto summary = seissol::statistics::parallelSummary(localPeak);
    logInfo() << "Peak host memory:" << UnitByte.formatPrefix(summary.mean).c_str()
              << "(per rank; range: [" << UnitByte.formatPrefix(summary.min).c_str() << ","
              << UnitByte.formatPrefix(summary.max).c_str() << "])";
    section["hostPeak"] = monitoring::toJson(summary);
  }
#ifdef ACL_DEVICE
  auto& api = device::DeviceInstance::instance().api();
  section["device"]["occupied"] = runSummary.overRanks(
      "deviceMemory" + moment, static_cast<double>(api.getCurrentlyOccupiedMem()));
  section["device"]["unified"] =
      monitoring::statisticsOverRanks(static_cast<double>(api.getCurrentlyOccupiedUnifiedMem()));
  section["device"]["capacity"] =
      monitoring::statisticsOverRanks(static_cast<double>(api.getMaxAvailableMem()));
#endif
  return section;
}

//! @brief The configurations the cells of the run compute in, by their names.
json::Value configurationNames(const seissol::SeisSol& seissolInstance) {
  auto names = json::Value::array();
  for (const auto config : seissolInstance.parameters().model.configs()) {
    names.push(configName(configValue(config)));
  }
  return names;
}

void closeSeisSol(seissol::SeisSol& seissolInstance) {
  logInfo() << "Cleaning up memory.";

  // deallocate memory manager
  seissolInstance.deleteMemoryManager();
}

} // namespace

void seissolMain(seissol::SeisSol& seissolInstance) {
  createOutputDirectory(seissolInstance);

  auto& summary = seissolInstance.runSummary();
  summary["configurations"] = configurationNames(seissolInstance);
  summary.start(seissolInstance.parameters().output.prefix, seissolInstance.runId());

  seissol::Stopwatch initWatch;
  initWatch.start();
  initSeisSol(seissolInstance);
  reportHardwareRelatedStatus(seissolInstance);

  // just put a barrier here to make sure everyone is synched
  logInfo() << "Finishing initialization...";
  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());
  summary.setRankValue("initTotal", initWatch.pause());
  summary["time"]["init"]["total"] =
      monitoring::toJson(initWatch.printTime("Initialization done in:"));
  summary["memory"]["afterInit"] = memorySection(summary, "AfterInit");

  seissol::Stopwatch watch;
  logInfo() << "Starting simulation.";
  summary.setStage("simulation");
  watch.start();
  seissolInstance.simulator().simulate(seissolInstance);
  watch.pause();
  static_cast<void>(watch.printTime("Time spent in simulation:"));

  // make sure everyone is really done
  logInfo() << "Simulation done.";
  seissol::Mpi::barrier(seissol::Mpi::mpi.comm());
  DeviceTimeline::stopKernelTrace();
  monitoring::Trace::write(seissolInstance.parameters().output.prefix + "-trace.json",
                           seissol::Mpi::mpi.comm());
  summary["memory"]["atEnd"] = memorySection(summary, "AtEnd");
  summary.writeRanks();
  summary.finish();

  closeSeisSol(seissolInstance);
}

} // namespace seissol::initializer::initprocedure

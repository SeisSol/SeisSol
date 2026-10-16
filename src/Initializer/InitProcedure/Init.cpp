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
#include "Monitoring/Memory.h"
#include "Monitoring/RunSummary.h"
#include "Monitoring/Stopwatch.h"
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
  const auto rank = seissol::Mpi::mpi.rank();
  std::size_t copyCells = 0;
  for (const auto& element : mesh.getElements()) {
    if (geometry::isCopy(element, rank)) {
      ++copyCells;
    }
  }
  json::Value section;
  section["cells"] =
      monitoring::statisticsOverRanks(static_cast<double>(mesh.getElements().size()));
  section["copyCells"] = monitoring::statisticsOverRanks(static_cast<double>(copyCells));
  section["ghostCells"] =
      monitoring::statisticsOverRanks(static_cast<double>(mesh.linearGhostlayer().size()));
  section["faultFaces"] =
      monitoring::statisticsOverRanks(static_cast<double>(mesh.getFault().size()));
  return section;
}

void initSeisSol(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.parameters();
  auto& summary = seissolInstance.runSummary();

  // set g
  seissolInstance.gravitationSetup().acceleration = seissolParams.model.gravitationalAcceleration;

  // the stages that do not time themselves are timed here; the summary gets all of them
  const auto timed = [&](const char* name, const char* message, auto&& stage) {
    seissol::Stopwatch watch;
    watch.start();
    stage(seissolInstance);
    watch.pause();
    summary["time"]["init"][name] = monitoring::toJson(watch.printTime(message));
  };

  // initialization procedure
  summary.setStage("init.mesh");
  seissol::initializer::initprocedure::initMesh(seissolInstance);
  summary.setStage("init.layout");
  timed("layout", "Layout initialized in:", seissol::initializer::initprocedure::initLayout);
  summary["mesh"] = meshSection(seissolInstance);
  summary.setStage("init.model");
  seissol::initializer::initprocedure::initModel(seissolInstance);
  summary.setStage("init.sideConditions");
  timed("sideConditions",
        "Side conditions initialized in:",
        seissol::initializer::initprocedure::initSideConditions);
  summary.setStage("init.io");
  timed("io", "Output initialized in:", seissol::initializer::initprocedure::initIO);

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
  parallel["asyncMode"] = [] {
    switch (async::Config::mode()) {
    case async::Mode::Thread:
      return "thread";
    case async::Mode::MPI:
      return "mpi";
    default:
      return "sync";
    }
  }();
}

/**
 * @brief The memory the run holds, per rank: the peak of the host memory, and what is allocated on
 * the device. Collective; prints the host peak as well.
 */
json::Value memorySection() {
  json::Value section;
  const auto peak = monitoring::peakHostMemory();
  if (peak.has_value()) {
    const auto summary = seissol::statistics::parallelSummary(static_cast<double>(peak.value()));
    logInfo() << "Peak host memory:" << UnitByte.formatPrefix(summary.mean).c_str()
              << "(per rank; range: [" << UnitByte.formatPrefix(summary.min).c_str() << ","
              << UnitByte.formatPrefix(summary.max).c_str() << "])";
    section["hostPeak"] = monitoring::toJson(summary);
  }
#ifdef ACL_DEVICE
  auto& api = device::DeviceInstance::instance().api();
  section["device"]["occupied"] =
      monitoring::statisticsOverRanks(static_cast<double>(api.getCurrentlyOccupiedMem()));
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
  initWatch.pause();
  summary["time"]["init"]["total"] =
      monitoring::toJson(initWatch.printTime("Initialization done in:"));
  summary["memory"]["afterInit"] = memorySection();

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
  summary["memory"]["atEnd"] = memorySection();
  summary.finish();

  closeSeisSol(seissolInstance);
}

} // namespace seissol::initializer::initprocedure

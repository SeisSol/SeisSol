// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Init.h"

#include <Memory/Tree/Layer.h>
#include <Monitoring/Stopwatch.h>
#include <utils/logger.h>

#include "InitIO.h"
#include "InitMesh.h"
#include "InitModel.h"
#include "InitSideConditions.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Parallel/MPI.h"
#include "ResultWriter/ThreadsPinningWriter.h"
#include "SeisSol.h"

#ifdef ACL_DEVICE
#include "Monitoring/Unit.h"
#include "Numerical/Statistics.h"
#include <ostream>
#include <sstream>
#endif

namespace {

void reportDeviceMemoryStatus() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  if (device.api->getCurrentlyOccupiedMem() > device.api->getMaxAvailableMem()) {
    std::stringstream stream;

    stream << "Memory of the device is overloaded." << std::endl
           << "Totally allocated device memory: "
           << UnitByte.formatPrefix(device.api->getCurrentlyOccupiedMem()) << std::endl
           << "Allocated unified memory: "
           << UnitByte.formatPrefix(device.api->getCurrentlyOccupiedUnifiedMem()) << std::endl
           << "Memory capacity of device: "
           << UnitByte.formatPrefix(device.api->getMaxAvailableMem());

    logError() << stream.str();
  } else {
    const double fraction = device.api->getCurrentlyOccupiedMem() /
                            static_cast<double>(device.api->getMaxAvailableMem());
    const auto summary = seissol::statistics::parallelSummary(fraction * 100.0);
    logInfo() << "occupied memory on devices (%):"
              << " mean =" << summary.mean << " std =" << summary.std << " min =" << summary.min
              << " median =" << summary.median << " max =" << summary.max;
  }
#endif
}

void initSeisSol(seissol::SeisSol& seissolInstance) {
  const auto& seissolParams = seissolInstance.getSeisSolParameters();

  // set g
  seissolInstance.getGravitationSetup().acceleration =
      seissolParams.model.gravitationalAcceleration;

  // initialization procedure
  seissol::initializer::initprocedure::initMesh(seissolInstance);
  seissol::initializer::initprocedure::initModel(seissolInstance);
  seissol::initializer::initprocedure::initSideConditions(seissolInstance);
  seissol::initializer::initprocedure::initIO(seissolInstance);

  // synchronize data to device
  seissolInstance.getMemoryManager().synchronizeTo(seissol::initializer::AllocationPlace::Device);

  // set up simulator
  auto& sim = seissolInstance.simulator();
  sim.setUsePlasticity(seissolParams.model.plasticity);
  sim.setFinalTime(seissolParams.timeStepping.endTime);
}

void reportHardwareRelatedStatus(seissol::SeisSol& seissolInstance) {
  reportDeviceMemoryStatus();

  const auto& seissolParams = seissolInstance.getSeisSolParameters();
  writer::ThreadsPinningWriter pinningWriter(seissolParams.output.prefix);
  pinningWriter.write(seissolInstance.getPinning());
}

void closeSeisSol(seissol::SeisSol& seissolInstance) {
  logInfo() << "Closing IO.";
  // cleanup IO
  seissolInstance.waveFieldWriter().close();
  seissolInstance.faultWriter().close();
  seissolInstance.freeSurfaceWriter().close();

  // deallocate memory manager
  seissolInstance.deleteMemoryManager();
}

} // namespace

void seissol::initializer::initprocedure::seissolMain(seissol::SeisSol& seissolInstance) {
  initSeisSol(seissolInstance);
  reportHardwareRelatedStatus(seissolInstance);

  // just put a barrier here to make sure everyone is synched
  logInfo() << "Finishing initialization...";
  seissol::MPI::barrier(seissol::MPI::mpi.comm());

  seissol::Stopwatch watch;
  logInfo() << "Starting simulation.";
  watch.start();
  seissolInstance.simulator().simulate(seissolInstance);
  watch.pause();
  watch.printTime("Time spent in simulation:");

  // make sure everyone is really done
  logInfo() << "Simulation done.";
  seissol::MPI::barrier(seissol::MPI::mpi.comm());

  closeSeisSol(seissolInstance);
}

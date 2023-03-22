#include "Init.hpp"
#include "InitCells.hpp"
#include "InitIO.hpp"
#include "InitLts.hpp"
#include "InitMesh.hpp"
#include "InitSideConditions.hpp"
#include "SeisSol.h"
#include "Initializer/InputParameters.hpp"

void reportDeviceMemoryStatus() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  constexpr size_t GB = 1024 * 1024 * 1024;
  const auto rank = seissol::MPI::mpi.rank();
  if (device.api->getCurrentlyOccupiedMem() > device.api->getMaxAvailableMem()) {
    std::stringstream stream;

    stream << "Device(" << rank << ")  memory is overloaded.\n"
           << "Totally allocated device memory, GB: " << device.api->getCurrentlyOccupiedMem() / GB << '\n'
           << "Allocated unified memory, GB: " << device.api->getCurrentlyOccupiedUnifiedMem() / GB << '\n'
           << "Memory capacity of device, GB: " << device.api->getMaxAvailableMem() / GB;

    logError() << stream.str();
  }
  else {
    double fraction = device.api->getCurrentlyOccupiedMem() / static_cast<double>(device.api->getMaxAvailableMem());
    logInfo() << "occupied memory on device(" << rank << "): " << fraction * 100.0 << "%";
  }
#endif
}

void initSeisSol() {
    const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();
    seissol::initializer::initprocedure::LtsInfo ltsInfo; // TODO: make ltsInfo obsolete

    // set constants

    // set g
    seissol::SeisSol::main.getGravitationSetup().acceleration = ssp.model.gravitationalAcceleration;

    // initialization procedure
    seissol::initializer::initprocedure::initIOPreLts();
    seissol::initializer::initprocedure::initMesh();
    seissol::initializer::initprocedure::initLts(ltsInfo);
    seissol::initializer::initprocedure::initCells(ltsInfo);
    seissol::initializer::initprocedure::initSideConditions(ltsInfo);
    seissol::initializer::initprocedure::initIOPostLts(ltsInfo);

    // set up simulator
    auto& sim = seissol::SeisSol::main.simulator();
    sim.setUsePlasticity(ssp.model.plasticity ? 1 : 0);
    sim.setFinalTime(ssp.end.endTime);

    // report status (TODO: move somewhere better)
    reportDeviceMemoryStatus();
}

void closeSeisSol() {
  logInfo() << "Closing IO.";
    // cleanup IO
	seissol::SeisSol::main.waveFieldWriter().close();
	seissol::SeisSol::main.checkPointManager().close();
	seissol::SeisSol::main.faultWriter().close();
	seissol::SeisSol::main.freeSurfaceWriter().close();

    // deallocate memory manager
    seissol::SeisSol::main.deleteMemoryManager();
}

void seissol::initializer::initprocedure::seissolMain() {
  initSeisSol();
  logInfo() << "Starting simulation.";
  seissol::SeisSol::main.simulator().simulate();
  logInfo() << "Simulation done.";
  closeSeisSol();
}

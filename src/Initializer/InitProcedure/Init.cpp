#include "Init.hpp"
#include "InitModel.hpp"
#include "InitIO.hpp"
#include "InitMesh.hpp"
#include "InitSideConditions.hpp"
#include "SeisSol.h"
#include "Initializer/InputParameters.hpp"
#include "Parallel/MPI.h"

static void reportDeviceMemoryStatus() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  constexpr size_t GB = 1024 * 1024 * 1024;
  const auto rank = seissol::MPI::mpi.rank();
  if (device.api->getCurrentlyOccupiedMem() > device.api->getMaxAvailableMem()) {
    std::stringstream stream;

    stream << "Device(" << rank << ")  memory is overloaded.\n"
           << "Totally allocated device memory, GB: " << device.api->getCurrentlyOccupiedMem() / GB
           << '\n'
           << "Allocated unified memory, GB: " << device.api->getCurrentlyOccupiedUnifiedMem() / GB
           << '\n'
           << "Memory capacity of device, GB: " << device.api->getMaxAvailableMem() / GB;

    logError() << stream.str();
  } else {
    double fraction = device.api->getCurrentlyOccupiedMem() /
                      static_cast<double>(device.api->getMaxAvailableMem());
    logInfo() << "occupied memory on device(" << rank << "): " << fraction * 100.0 << "%";
  }
#endif
}

static void initSeisSol() {
  const auto& ssp = seissol::SeisSol::main.getSeisSolParameters();

  // set g
  seissol::SeisSol::main.getGravitationSetup().acceleration = ssp.model.gravitationalAcceleration;

  // initialization procedure
  seissol::initializer::initprocedure::initMesh();
  seissol::initializer::initprocedure::initModel();
  seissol::initializer::initprocedure::initSideConditions();
  seissol::initializer::initprocedure::initIO();

  // set up simulator
  auto& sim = seissol::SeisSol::main.simulator();
  sim.setUsePlasticity(ssp.model.plasticity);
  sim.setFinalTime(ssp.end.endTime);

  // report status (TODO(David): move somewhere better)
  reportDeviceMemoryStatus();
}

static void closeSeisSol() {
  logInfo(seissol::MPI::mpi.rank()) << "Closing IO.";
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

  // just put a barrier here to make sure everyone is synched
  logInfo(seissol::MPI::mpi.rank()) << "Finishing initialization...";
  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());

  seissol::Stopwatch watch;
  logInfo(seissol::MPI::mpi.rank()) << "Starting simulation.";
  watch.start();
  seissol::SeisSol::main.simulator().simulate();
  watch.pause();
  watch.printTime("Time spent in simulation:");

  // make sure everyone is really done
  logInfo(seissol::MPI::mpi.rank()) << "Simulation done.";
  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());

  closeSeisSol();
}

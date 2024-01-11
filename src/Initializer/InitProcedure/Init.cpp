#include "Init.hpp"
#include "InitModel.hpp"
#include "InitIO.hpp"
#include "InitMesh.hpp"
#include "InitSideConditions.hpp"
#include "SeisSol.h"
#include "Initializer/InputParameters.hpp"
#include "Parallel/MPI.h"
#include <Numerical_aux/Statistics.h>
#include "ResultWriter/ThreadsPinningWriter.h"
#include <sstream>
#include "Monitoring/Unit.hpp"

namespace {

static void reportDeviceMemoryStatus() {
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
  const auto rank = seissol::MPI::mpi.rank();
  if (device.api->getCurrentlyOccupiedMem() > device.api->getMaxAvailableMem()) {
    std::stringstream stream;

    stream << "Memory of device (" << rank << ") is overloaded." << std::endl
           << "Totally allocated device memory: "
           << UnitByte.formatPrefix(device.api->getCurrentlyOccupiedMem()) << std::endl
           << "Allocated unified memory: "
           << UnitByte.formatPrefix(device.api->getCurrentlyOccupiedUnifiedMem()) << std::endl
           << "Memory capacity of device: "
           << UnitByte.formatPrefix(device.api->getMaxAvailableMem());

    logError() << stream.str();
  } else {
    double fraction = device.api->getCurrentlyOccupiedMem() /
                      static_cast<double>(device.api->getMaxAvailableMem());
    const auto summary = seissol::statistics::parallelSummary(fraction * 100.0);
    logInfo(rank) << "occupied memory on devices (%):"
                  << " mean =" << summary.mean << " std =" << summary.std << " min =" << summary.min
                  << " median =" << summary.median << " max =" << summary.max;
  }
#endif
}

static void initSeisSol() {
  const auto& seissolParams = seissol::SeisSol::main.getSeisSolParameters();

  // set g
  seissol::SeisSol::main.getGravitationSetup().acceleration =
      seissolParams.model.gravitationalAcceleration;

  // initialization procedure
  seissol::initializer::initprocedure::initMesh();
  seissol::initializer::initprocedure::initModel();
  seissol::initializer::initprocedure::initSideConditions();
  seissol::initializer::initprocedure::initIO();

  // set up simulator
  auto& sim = seissol::SeisSol::main.simulator();
  sim.setUsePlasticity(seissolParams.model.plasticity);
  sim.setFinalTime(seissolParams.end.endTime);
}

static void reportHardwareRelatedStatus() {
  reportDeviceMemoryStatus();

  const auto& seissolParams = SeisSol::main.getSeisSolParameters();
  writer::ThreadsPinningWriter pinningWriter(seissolParams.output.prefix);
  pinningWriter.write(SeisSol::main.getPinning());
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

} // namespace

void seissol::initializer::initprocedure::seissolMain() {
  initSeisSol();
  reportHardwareRelatedStatus();

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

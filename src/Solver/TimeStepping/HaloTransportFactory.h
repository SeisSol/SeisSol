// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORTFACTORY_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORTFACTORY_H_

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <memory>

namespace seissol::solver {

class TaskHaloEngine;

/**
 * Creates the transports of one process, and holds what they share.
 */
class HaloTransportFactory {
  public:
  /**
   * Persistent MPI transports set up their requests once and restart them for each transfer.
   */
  HaloTransportFactory(Mpi::DataTransferMode mode, bool persistent);
  ~HaloTransportFactory();

  HaloTransportFactory(const HaloTransportFactory&) = delete;
  HaloTransportFactory(HaloTransportFactory&&) = delete;
  HaloTransportFactory& operator=(const HaloTransportFactory&) = delete;
  HaloTransportFactory& operator=(HaloTransportFactory&&) = delete;

  /**
   * Creates the transport of the regions between a copy layer and a ghost layer. The transports
   * need to be destroyed before the factory.
   */
  std::unique_ptr<HaloTransport> create(const RemoteClusterPair& regions);

  /**
   * What starts and tests the exchanges, with a task graph on the host (see TaskHaloTransport);
   * null otherwise.
   */
  [[nodiscard]] TaskHaloEngine* taskEngine() { return taskEngine_.get(); }

  private:
  Mpi::DataTransferMode mode_;
  bool persistent_;
  std::unique_ptr<TaskHaloEngine> taskEngine_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORTFACTORY_H_

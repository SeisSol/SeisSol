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

#include <cstddef>
#include <memory>

namespace seissol::time_stepping {

class ExchangeScheduler;

/**
 * Creates the transports of one process, and holds what they share.
 */
class HaloTransportFactory {
  public:
  /**
   * Persistent MPI transports set up their requests once and restart them for each transfer. The
   * transports on device streams (CCL, stream-aware MPI, SHMEM) set up their library here,
   * collectively over all processes: for all exchanges in an order that is the same on all
   * processes, or with `perDirection`, for each direction between two time clusters.
   */
  HaloTransportFactory(Mpi::DataTransferMode mode,
                       bool persistent,
                       bool perDirection,
                       std::size_t clusterCount);
  ~HaloTransportFactory();

  HaloTransportFactory(const HaloTransportFactory&) = delete;
  HaloTransportFactory(HaloTransportFactory&&) = delete;
  HaloTransportFactory& operator=(const HaloTransportFactory&) = delete;
  HaloTransportFactory& operator=(HaloTransportFactory&&) = delete;

  /**
   * Creates the transport between the copy layer of `cluster` and the ghost layer of
   * `otherCluster`. The transports need to be destroyed before the factory.
   */
  std::unique_ptr<HaloTransport> create(const solver::RemoteClusterPair& regions,
                                        std::size_t cluster,
                                        std::size_t otherCluster);

  /**
   * The scheduler of the CCL transports; null for the others.
   */
  [[nodiscard]] ExchangeScheduler* scheduler() { return scheduler_.get(); }

  private:
  Mpi::DataTransferMode mode_;
  bool persistent_;
  std::unique_ptr<ExchangeScheduler> scheduler_;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORTFACTORY_H_

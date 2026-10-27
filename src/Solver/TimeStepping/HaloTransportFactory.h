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
#include <vector>

namespace seissol::solver {

class ExchangeScheduler;
class TaskHaloEngine;

/**
 * Creates the transports of one process, and holds what they share.
 */
class HaloTransportFactory {
  public:
  /**
   * Persistent MPI transports set up their requests once and restart them for each transfer. The
   * transports on device streams (CCL, stream-aware MPI, SHMEM) set up their library here,
   * collectively over all processes: for all exchanges in an order that is the same on all
   * processes (on `lanes` lanes, see StreamExchangeScheduler), or with `perDirection`, for each
   * direction between two clusters. The clusters are the time clusters of each of `configCount`
   * configurations (see ExchangeScheduler).
   */
  HaloTransportFactory(Mpi::DataTransferMode mode,
                       bool persistent,
                       bool perDirection,
                       std::size_t clusterCount,
                       std::size_t configCount,
                       std::size_t lanes = 1);
  ~HaloTransportFactory();

  HaloTransportFactory(const HaloTransportFactory&) = delete;
  HaloTransportFactory(HaloTransportFactory&&) = delete;
  HaloTransportFactory& operator=(const HaloTransportFactory&) = delete;
  HaloTransportFactory& operator=(HaloTransportFactory&&) = delete;

  /**
   * Creates the transport between the copy layer of `cluster` and the ghost layer of
   * `otherCluster`. The transports need to be destroyed before the factory.
   */
  std::unique_ptr<HaloTransport>
      create(const RemoteClusterPair& regions, std::size_t cluster, std::size_t otherCluster);

  /**
   * The scheduler of the CCL transports; null for the others.
   */
  [[nodiscard]] ExchangeScheduler* scheduler() { return scheduler_.get(); }

  /**
   * Lets the transports created from now on that run on device streams (ccl, stream-mpi, shmem)
   * move only the regions with the processes for which `onNode` is true (by their rank), and the
   * others with GPU-aware MPI from the host (see SplitHaloTransport); experimental.
   */
  void splitByNode(std::vector<bool> onNode);

  /// whether some transports move regions with MPI from the host (see splitByNode())
  [[nodiscard]] bool split() const { return !onNode_.empty(); }

  /**
   * What starts and tests the exchanges, with a task graph on the host (see TaskHaloTransport);
   * null otherwise.
   */
  [[nodiscard]] TaskHaloEngine* taskEngine() { return taskEngine_.get(); }

  private:
  Mpi::DataTransferMode mode_;
  bool persistent_;
  // by rank, whether the regions with the process stay on the device streams; empty for all
  std::vector<bool> onNode_;
  std::unique_ptr<ExchangeScheduler> scheduler_;
  std::unique_ptr<TaskHaloEngine> taskEngine_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOTRANSPORTFACTORY_H_

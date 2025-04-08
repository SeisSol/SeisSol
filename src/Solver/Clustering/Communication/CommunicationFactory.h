// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_COMMUNICATIONFACTORY_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_COMMUNICATIONFACTORY_H_

#include "Parallel/MPI.h"
#include "memory"
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/Communication/CCLNeighborCluster.h>
#include <Solver/Clustering/Communication/DirectMPINeighborCluster.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>

namespace seissol::solver::clustering::communication {
enum class CommunicationMode {
  DirectMPI,
  HostSurrogateMPI,
  DirectCCL,
};

class CommunicationClusterFactory {
  private:
  std::vector<void*> comms;
  std::size_t globalClusterCount;
  CommunicationMode mode;

  public:
  CommunicationClusterFactory(CommunicationMode mode, std::size_t clusterCount);

  void prepare();

  std::shared_ptr<SendNeighborCluster>
      getSend(std::size_t cluster,
              const std::vector<RemoteCluster>& remoteClusters,
              const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
              double priority);

  std::shared_ptr<RecvNeighborCluster>
      getRecv(std::size_t cluster,
              const std::vector<RemoteCluster>& remoteClusters,
              const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
              double priority);

  std::vector<std::shared_ptr<RecvNeighborCluster>>
      getAllRecvs(const HaloCommunication& comm,
                  const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                  double priority);

  std::vector<std::shared_ptr<SendNeighborCluster>>
      getAllSends(const HaloCommunication& comm,
                  const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                  double priority);
};
} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_COMMUNICATIONFACTORY_H_

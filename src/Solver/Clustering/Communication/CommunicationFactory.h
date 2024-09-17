// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Parallel/MPI.h"
#include "memory"
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/Communication/CCLNeighborCluster.h>
#include <Solver/Clustering/Communication/DirectMPINeighborCluster.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>

#if defined(ACL_DEVICE) && defined(USE_CCL)
#include "CCLCluster.h"
#include "CCLSetup.h"
#endif

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

  std::unique_ptr<SendNeighborCluster> getSend(std::size_t cluster,
                                               const std::vector<RemoteCluster>& remoteClusters);

  std::unique_ptr<RecvNeighborCluster> getRecv(std::size_t cluster,
                                               const std::vector<RemoteCluster>& remoteClusters);

  std::vector<std::unique_ptr<RecvNeighborCluster>> getAllRecvs(const HaloCommunication& comm);

  std::vector<std::unique_ptr<RecvNeighborCluster>> getAllSends(const HaloCommunication& comm);
};
} // namespace seissol::solver::clustering::communication

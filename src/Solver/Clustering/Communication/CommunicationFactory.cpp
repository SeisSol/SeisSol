// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "CommunicationFactory.h"

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

CommunicationClusterFactory::CommunicationClusterFactory(CommunicationMode mode,
                                                         std::size_t clusterCount)
    : mode(mode), globalClusterCount(clusterCount) {}

void CommunicationClusterFactory::prepare() {
#ifdef ACL_DEVICE
  if (mode == CommunicationMode::DirectCCL) {
    comms = seissol::solver::clustering::communication::createComms(globalClusterCount);
  }
#endif
}

std::unique_ptr<SendNeighborCluster>
    CommunicationClusterFactory::getSend(std::size_t cluster,
                                         const std::vector<RemoteCluster>& remoteClusters) {
  switch (mode) {
#ifdef ACL_DEVICE
  case CommunicationMode::HostSurrogateMPI: {
    // return std::make_unique<>(remoteClusters);
    return nullptr;
  }
#ifdef USE_CCL
  case CommunicationMode::DirectCCL: {
    return std::make_unique<CCLSendNeighborCluster>(remoteClusters, comms.at(cluster));
  }
#endif
#endif // ACL_DEVICE
  case CommunicationMode::DirectMPI: {
    return std::make_unique<DirectMPISendNeighborCluster>(remoteClusters);
  }
  default: {
    return nullptr;
  }
  }
}

std::unique_ptr<RecvNeighborCluster>
    CommunicationClusterFactory::getRecv(std::size_t cluster,
                                         const std::vector<RemoteCluster>& remoteClusters) {
  switch (mode) {
#ifdef ACL_DEVICE
  case CommunicationMode::HostSurrogateMPI: {
    // return std::make_unique<>(remoteClusters);
    return nullptr;
  }
#ifdef USE_CCL
  case CommunicationMode::DirectCCL: {
    return std::make_unique<CCLRecvNeighborCluster>(remoteClusters, comms.at(cluster));
  }
#endif
#endif // ACL_DEVICE
  case CommunicationMode::DirectMPI: {
    return std::make_unique<DirectMPIRecvNeighborCluster>(remoteClusters);
  }
  default: {
    return nullptr;
  }
  }
}

std::vector<std::unique_ptr<RecvNeighborCluster>>
    CommunicationClusterFactory::getAllRecvs(const HaloCommunication& comm) {
  std::vector<std::unique_ptr<RecvNeighborCluster>> clusters;
  clusters.reserve(comm.ghost.size());
  for (std::size_t i = 0; i < comm.ghost.size(); ++i) {
    clusters.emplace_back(getRecv(i, comm.ghost.at(i)));
  }
  return clusters;
}

std::vector<std::unique_ptr<RecvNeighborCluster>>
    CommunicationClusterFactory::getAllSends(const HaloCommunication& comm) {
  std::vector<std::unique_ptr<RecvNeighborCluster>> clusters;
  clusters.reserve(comm.copy.size());
  for (std::size_t i = 0; i < comm.ghost.size(); ++i) {
    clusters.emplace_back(getSend(i, comm.copy.at(i)));
  }
  return clusters;
}
} // namespace seissol::solver::clustering::communication

// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CommunicationFactory.h"

#include "Parallel/MPI.h"
#include "memory"
#include <Parallel/Host/CpuExecutor.h>
#include <Solver/Clustering/AbstractTimeCluster.h>
#include <Solver/Clustering/Communication/CCLNeighborCluster.h>
#include <Solver/Clustering/Communication/DirectMPINeighborCluster.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <memory>

#if defined(ACL_DEVICE) && defined(USE_CCL)
#include "CCLNeighborCluster.h"
#include "CCLSetup.h"
#endif

namespace seissol::solver::clustering::communication {

CommunicationClusterFactory::CommunicationClusterFactory(CommunicationMode mode,
                                                         std::size_t clusterCount)
    : mode(mode), globalClusterCount(clusterCount) {}

void CommunicationClusterFactory::prepare() {
#if defined(ACL_DEVICE) && defined(USE_CCL)
  if (mode == CommunicationMode::DirectCCL) {
    comms = seissol::solver::clustering::communication::createComms(globalClusterCount);
  }
#endif
}

std::shared_ptr<SendNeighborCluster> CommunicationClusterFactory::getSend(
    std::size_t cluster,
    const std::vector<RemoteCluster>& remoteClusters,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority) {
  switch (mode) {
#ifdef ACL_DEVICE
  case CommunicationMode::HostSurrogateMPI: {
    // return std::make_shared<>(remoteClusters);
    return nullptr;
  }
#ifdef USE_CCL
  case CommunicationMode::DirectCCL: {
    return std::make_shared<CCLSendNeighborCluster>(
        remoteClusters, comms.at(cluster), cpuExecutor, priority);
  }
#endif
#endif // ACL_DEVICE
  case CommunicationMode::DirectMPI: {
    return std::make_shared<DirectMPISendNeighborCluster>(remoteClusters, cpuExecutor, priority);
  }
  default: {
    return nullptr;
  }
  }
}

std::shared_ptr<RecvNeighborCluster> CommunicationClusterFactory::getRecv(
    std::size_t cluster,
    const std::vector<RemoteCluster>& remoteClusters,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority) {
  switch (mode) {
#ifdef ACL_DEVICE
  case CommunicationMode::HostSurrogateMPI: {
    // return std::make_shared<>(remoteClusters);
    return nullptr;
  }
#ifdef USE_CCL
  case CommunicationMode::DirectCCL: {
    return std::make_shared<CCLRecvNeighborCluster>(
        remoteClusters, comms.at(cluster), cpuExecutor, priority);
  }
#endif
#endif // ACL_DEVICE
  case CommunicationMode::DirectMPI: {
    return std::make_shared<DirectMPIRecvNeighborCluster>(remoteClusters, cpuExecutor, priority);
  }
  default: {
    return nullptr;
  }
  }
}

std::vector<std::shared_ptr<RecvNeighborCluster>> CommunicationClusterFactory::getAllRecvs(
    const HaloCommunication& comm,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority) {
  std::vector<std::shared_ptr<RecvNeighborCluster>> clusters;
  clusters.reserve(comm.ghost.size());
  for (std::size_t i = 0; i < comm.ghost.size(); ++i) {
    clusters.emplace_back(getRecv(i, comm.ghost.at(i), cpuExecutor, priority));
  }
  return clusters;
}

std::vector<std::shared_ptr<SendNeighborCluster>> CommunicationClusterFactory::getAllSends(
    const HaloCommunication& comm,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority) {
  std::vector<std::shared_ptr<SendNeighborCluster>> clusters;
  clusters.reserve(comm.copy.size());
  for (std::size_t i = 0; i < comm.ghost.size(); ++i) {
    clusters.emplace_back(getSend(i, comm.copy.at(i), cpuExecutor, priority));
  }
  return clusters;
}
} // namespace seissol::solver::clustering::communication

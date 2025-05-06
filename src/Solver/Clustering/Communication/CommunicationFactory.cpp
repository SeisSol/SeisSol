// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CommunicationFactory.h"

#include "Parallel/MPI.h"
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

std::shared_ptr<NeighborCluster> CommunicationClusterFactory::getPair(
    std::size_t cluster,
    double stepWidth,
    std::size_t stepRate,
    const std::vector<RemoteCluster>& remoteSend,
    const std::vector<RemoteCluster>& remoteRecv,
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
    return std::make_shared<CCLNeighborCluster>(
        stepWidth, stepRate, remoteSend, remoteRecv, comms.at(cluster), cpuExecutor, priority);
  }
#endif
#endif // ACL_DEVICE
  case CommunicationMode::DirectMPI: {
    return std::make_shared<DirectMPINeighborCluster>(
        stepWidth, stepRate, remoteSend, remoteRecv, cpuExecutor, priority);
  }
  default: {
    return nullptr;
  }
  }
}

std::vector<std::shared_ptr<NeighborCluster>> CommunicationClusterFactory::get(
    const HaloCommunication& comm,
    const initializer::ClusterLayout& layout,
    const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
    double priority) {
  std::vector<std::shared_ptr<NeighborCluster>> clusters;
  clusters.resize(comm.copy.size());

  std::size_t commCluster = 0;
  for (std::size_t i = 0; i < comm.ghost.size(); ++i) {
    clusters[i] = getPair(i,
                          layout.timestepRate(i),
                          layout.clusterRate(i),
                          comm.copy[i],
                          comm.ghost[i],
                          cpuExecutor,
                          priority);
  }
  return clusters;
}
} // namespace seissol::solver::clustering::communication

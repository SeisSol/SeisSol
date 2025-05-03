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
#include <Solver/Clustering/Communication/CCLNeighborCluster2.h>
#include <Solver/Clustering/Communication/DirectMPINeighborClusterBlocking.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <memory>

#if defined(ACL_DEVICE) && defined(USE_CCL)
#include "CCLNeighborCluster2.h"
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

std::pair<std::shared_ptr<SendNeighborCluster>, std::shared_ptr<RecvNeighborCluster>>
    CommunicationClusterFactory::getPair(
        std::size_t cluster,
        const std::vector<RemoteCluster>& remoteSend,
        const std::vector<RemoteCluster>& remoteRecv,
        const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
        double priority) {
  switch (mode) {
#ifdef ACL_DEVICE
  case CommunicationMode::HostSurrogateMPI: {
    // return std::make_shared<>(remoteClusters);
    return {nullptr, nullptr};
  }
#ifdef USE_CCL
  case CommunicationMode::DirectCCL: {
    return {std::make_shared<CCLNeighborCluster>(
                remoteSend, remoteRecv, comms.at(cluster), cpuExecutor, priority),
            std::make_shared<CCLNeighborDummyCluster>(cpuExecutor, priority)};
  }
#endif
#endif // ACL_DEVICE
  case CommunicationMode::DirectMPI: {
    return {
        std::make_shared<DirectMPISendNeighborClusterBlocking>(remoteSend, cpuExecutor, priority),
        std::make_shared<DirectMPIRecvNeighborClusterBlocking>(remoteRecv, cpuExecutor, priority)};
  }
  default: {
    return {nullptr, nullptr};
  }
  }
}

std::pair<std::vector<std::shared_ptr<SendNeighborCluster>>,
          std::vector<std::shared_ptr<RecvNeighborCluster>>>
    CommunicationClusterFactory::get(
        const HaloCommunication& comm,
        const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
        double priority) {
  std::vector<std::shared_ptr<RecvNeighborCluster>> clustersRecv;
  std::vector<std::shared_ptr<SendNeighborCluster>> clustersSend;
  clustersSend.reserve(comm.copy.size());
  clustersRecv.reserve(comm.ghost.size());
  for (std::size_t i = 0; i < comm.ghost.size(); ++i) {
    const auto sendrecv = getPair(i, comm.copy.at(i), comm.ghost.at(i), cpuExecutor, priority);
    clustersSend.emplace_back(sendrecv.first);
    clustersRecv.emplace_back(sendrecv.second);
  }
  return {clustersSend, clustersRecv};
}
} // namespace seissol::solver::clustering::communication

// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CommunicationManager.h"

#include "Parallel/Pin.h"
#include <memory>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include "device.h"
#endif // ACL_DEVICE

namespace seissol::time_stepping {

AbstractCommunicationManager::AbstractCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters)
    : ghostClusters(std::move(ghostClusters)) {}
void AbstractCommunicationManager::reset(double newSyncTime) {
  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->setSyncTime(newSyncTime);

    // dereference first due to clang-tidy recommendation
    (*ghostCluster).reset();
  }
}

std::vector<std::unique_ptr<AbstractGhostTimeCluster>>*
    AbstractCommunicationManager::getGhostClusters() {
  return &ghostClusters;
}

bool AbstractCommunicationManager::poll() {
  bool finished = true;
  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->act();
    finished = finished && ghostCluster->synced();
  }
  return finished;
}

SerialCommunicationManager::SerialCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters)
    : AbstractCommunicationManager(std::move(ghostClusters)) {}

bool SerialCommunicationManager::checkIfFinished() const {
  for (const auto& ghostCluster : ghostClusters) {
    if (!ghostCluster->synced()) {
      return false;
    }
  }
  return true;
}

void SerialCommunicationManager::progression() { poll(); }

ThreadedCommunicationManager::ThreadedCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters,
    const seissol::parallel::Pinning* pinning)
    : AbstractCommunicationManager(std::move(ghostClusters)),
      helper([&]() { return this->poll(); }, pinning) {}

void ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool ThreadedCommunicationManager::checkIfFinished() const { return helper.finished(); }

void ThreadedCommunicationManager::reset(double newSyncTime) {
  helper.stop();
  AbstractCommunicationManager::reset(newSyncTime);
  helper.start();
}

} // namespace seissol::time_stepping

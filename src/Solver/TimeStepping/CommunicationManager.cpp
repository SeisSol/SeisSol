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
    ghostCluster->reset();
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
    : AbstractCommunicationManager(std::move(ghostClusters)), shouldReset(false), isFinished(false),
      pinning(pinning) {}

void ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool ThreadedCommunicationManager::checkIfFinished() const { return isFinished.load(); }

void ThreadedCommunicationManager::reset(double newSyncTime) {
  // Send signal to comm. thread to finish and wait.
  shouldReset.store(true);
  if (thread.joinable()) {
    thread.join();
  }

  // Reset flags and reset ghost clusters
  shouldReset.store(false);
  isFinished.store(false);
  AbstractCommunicationManager::reset(newSyncTime);

  // Start a new communication thread.
  // Note: Easier than keeping one alive, and not that expensive.
  thread = std::thread([this]() {
#ifdef ACL_DEVICE
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    device.api->setDevice(0);
#endif // ACL_DEVICE
    // Pin this thread to the last core
    // We compute the mask outside the thread because otherwise
    // it confuses profilers and debuggers!
    pinning->pinToFreeCPUs();
    while (!shouldReset.load() && !isFinished.load()) {
      isFinished.store(this->poll());
    }
  });
}

ThreadedCommunicationManager::~ThreadedCommunicationManager() {
  if (thread.joinable()) {
    thread.join();
  }
}

} // namespace seissol::time_stepping

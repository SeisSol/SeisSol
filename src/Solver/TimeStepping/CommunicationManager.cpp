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

namespace seissol::time_stepping {

AbstractCommunicationManager::AbstractCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters)
    : ghostClusters_(std::move(ghostClusters)) {}
void AbstractCommunicationManager::reset(double newSyncTime) {
  for (auto& ghostCluster : ghostClusters_) {
    ghostCluster->setSyncTime(newSyncTime);

    // dereference first due to clang-tidy recommendation
    (*ghostCluster).reset();
  }
}

std::vector<std::unique_ptr<AbstractGhostTimeCluster>>*
    AbstractCommunicationManager::getGhostClusters() {
  return &ghostClusters_;
}

bool AbstractCommunicationManager::poll() {
  bool finished = true;
  for (auto& ghostCluster : ghostClusters_) {
    ghostCluster->act();
    finished = finished && ghostCluster->synced();
  }
  return finished;
}

SerialCommunicationManager::SerialCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters)
    : AbstractCommunicationManager(std::move(ghostClusters)) {}

bool SerialCommunicationManager::checkIfFinished() const {
  for (const auto& ghostCluster : ghostClusters_) {
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
      helper_([&]() { return this->poll(); }, pinning) {}

void ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool ThreadedCommunicationManager::checkIfFinished() const { return helper_.finished(); }

void ThreadedCommunicationManager::reset(double newSyncTime) {
  helper_.stop();
  AbstractCommunicationManager::reset(newSyncTime);
  helper_.start();
}

} // namespace seissol::time_stepping

#include "CommunicationManager.h"

#include "Parallel/Pin.h"
#include <memory>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include "device.h"
#endif // ACL_DEVICE

seissol::time_stepping::AbstractCommunicationManager::AbstractCommunicationManager(
    seissol::time_stepping::AbstractCommunicationManager::ghostClusters_t ghostClusters)
    : ghostClusters(std::move(ghostClusters)) {}
void seissol::time_stepping::AbstractCommunicationManager::reset(double newSyncTime) {
  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->setSyncTime(newSyncTime);
    ghostCluster->reset();
  }
}

std::vector<std::unique_ptr<seissol::time_stepping::AbstractGhostTimeCluster>>*
    seissol::time_stepping::AbstractCommunicationManager::getGhostClusters() {
  return &ghostClusters;
}

bool seissol::time_stepping::AbstractCommunicationManager::poll() {
  bool finished = true;
  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->act();
    finished = finished && ghostCluster->synchronized();
  }
  return finished;
}

seissol::time_stepping::SerialCommunicationManager::SerialCommunicationManager(
    seissol::time_stepping::AbstractCommunicationManager::ghostClusters_t ghostClusters)
    : AbstractCommunicationManager(std::move(ghostClusters)) {}

bool seissol::time_stepping::SerialCommunicationManager::checkIfFinished() const {
  for (auto& ghostCluster : ghostClusters) {
    if (!ghostCluster->synchronized())
      return false;
  }
  return true;
}

void seissol::time_stepping::SerialCommunicationManager::progression() { poll(); }

seissol::time_stepping::ThreadedCommunicationManager::ThreadedCommunicationManager(
    seissol::time_stepping::AbstractCommunicationManager::ghostClusters_t ghostClusters,
    const seissol::parallel::Pinning* pinning)
    : AbstractCommunicationManager(std::move(ghostClusters)),
      helper([this]() { return this->poll(); }, pinning) {}

void seissol::time_stepping::ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool seissol::time_stepping::ThreadedCommunicationManager::checkIfFinished() const {
  return helper.finished();
}

void seissol::time_stepping::ThreadedCommunicationManager::reset(double newSyncTime) {
  AbstractCommunicationManager::reset(newSyncTime);

  helper.restart();
}

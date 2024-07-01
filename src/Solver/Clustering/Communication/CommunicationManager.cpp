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
    AbstractCommunicationManager::ghostClusters_t ghostClusters)
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
    finished = finished && ghostCluster->synchronized();
  }
  return finished;
}

SerialCommunicationManager::SerialCommunicationManager(
    AbstractCommunicationManager::ghostClusters_t ghostClusters)
    : AbstractCommunicationManager(std::move(ghostClusters)) {}

bool SerialCommunicationManager::checkIfFinished() const {
  for (auto& ghostCluster : ghostClusters) {
    if (!ghostCluster->synchronized())
      return false;
  }
  return true;
}

void SerialCommunicationManager::progression() { poll(); }

ThreadedCommunicationManager::ThreadedCommunicationManager(
    AbstractCommunicationManager::ghostClusters_t ghostClusters,
    const seissol::parallel::Pinning* pinning)
    : AbstractCommunicationManager(std::move(ghostClusters)),
      helper([this]() { return this->poll(); }, pinning) {}

void ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool ThreadedCommunicationManager::checkIfFinished() const { return helper.finished(); }

void ThreadedCommunicationManager::reset(double newSyncTime) {
  AbstractCommunicationManager::reset(newSyncTime);

  helper.restart();
}

} // namespace seissol::time_stepping

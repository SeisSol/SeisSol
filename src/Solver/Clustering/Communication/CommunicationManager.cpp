#include "CommunicationManager.h"

#include "Parallel/Pin.h"
#include <memory>
#include <utility>
#include <vector>

#ifdef ACL_DEVICE
#include "device.h"
#endif // ACL_DEVICE

namespace seissol::solver::clustering::communication {

AbstractCommunicationManager::AbstractCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters)
    : ghostClusters(std::move(ghostClusters)) {}
void AbstractCommunicationManager::reset() {
}

AbstractCommunicationManager::GhostClustersT&
    AbstractCommunicationManager::getGhostClusters() {
  return ghostClusters;
}

bool AbstractCommunicationManager::poll() {
  bool finished = true;
  for (auto& ghostCluster : ghostClusters) {
    const auto clusterDone = ghostCluster->poll();
    finished = finished && clusterDone;
  }
  return finished;
}

SerialCommunicationManager::SerialCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters)
    : AbstractCommunicationManager(std::move(ghostClusters)) {}

bool SerialCommunicationManager::checkIfFinished() const {
  for (auto& ghostCluster : ghostClusters) {
    if (!ghostCluster->poll()) {
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
      helper([this]() { return this->poll(); }, pinning) {}

void ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool ThreadedCommunicationManager::checkIfFinished() const { return helper.finished(); }

void ThreadedCommunicationManager::reset() {
  AbstractCommunicationManager::reset();

  helper.restart();
}

} // namespace seissol::solver::clustering::communication

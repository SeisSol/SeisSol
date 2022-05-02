#include "CommunicationManager.h"

#include "Parallel/Pin.h"

seissol::time_stepping::AbstractCommunicationManager::AbstractCommunicationManager(
    seissol::time_stepping::AbstractCommunicationManager::ghostClusters_t ghostClusters) : ghostClusters(std::move(ghostClusters)) {

}
void seissol::time_stepping::AbstractCommunicationManager::reset(double newSyncTime) {
  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->setSyncTime(newSyncTime);
    ghostCluster->reset();
  }
}

bool seissol::time_stepping::AbstractCommunicationManager::poll() {
  bool finished = true;
  for (auto& ghostCluster : ghostClusters) {
    ghostCluster->act();
    finished = finished && ghostCluster->synced();
  }
  return finished;
}

seissol::time_stepping::SerialCommunicationManager::SerialCommunicationManager(
    seissol::time_stepping::AbstractCommunicationManager::ghostClusters_t ghostClusters)
    : AbstractCommunicationManager(std::move(ghostClusters)) {

}

bool seissol::time_stepping::SerialCommunicationManager::checkIfFinished() const {
  for (auto& ghostCluster : ghostClusters) {
    if (!ghostCluster->synced()) return false;
  }
  return true;
}

void seissol::time_stepping::SerialCommunicationManager::progression() {
  poll();
}

seissol::time_stepping::ThreadedCommunicationManager::ThreadedCommunicationManager(
    seissol::time_stepping::AbstractCommunicationManager::ghostClusters_t ghostClusters,
    const seissol::parallel::Pinning* pinning)
    : AbstractCommunicationManager(std::move(ghostClusters)),
      thread(),
      shouldReset(false),
      isFinished(false),
      pinning(pinning) {
}

void seissol::time_stepping::ThreadedCommunicationManager::progression() {
  // Do nothing: Thread takes care of that.
}

bool seissol::time_stepping::ThreadedCommunicationManager::checkIfFinished() const {
  return isFinished.load();
}

void seissol::time_stepping::ThreadedCommunicationManager::reset(double newSyncTime) {
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
  thread = std::thread([this](){
    // Pin this thread to the last core
    // We compute the mask outside the thread because otherwise
    // it confuses profilers and debuggers!
    pinning->pinToFreeCPUs();
    while(!shouldReset.load() && !isFinished.load()) {
      isFinished.store(this->poll());
    }
  });
}

seissol::time_stepping::ThreadedCommunicationManager::~ThreadedCommunicationManager() {
  if (thread.joinable()) {
    thread.join();
  }
}




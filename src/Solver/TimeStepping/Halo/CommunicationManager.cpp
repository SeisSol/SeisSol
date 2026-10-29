// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "CommunicationManager.h"

#include "Parallel/Pin.h"
#include "Parallel/Runtime/HostRuntime.h"
#include "Solver/TimeStepping/Halo/TaskHaloTransport.h"

#include <algorithm>
#include <atomic>
#include <memory>
#include <utility>
#include <vector>

namespace seissol::solver {

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

std::vector<std::unique_ptr<GhostCluster>>* AbstractCommunicationManager::getGhostClusters() {
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

TaskCommunicationManager::TaskCommunicationManager(
    AbstractCommunicationManager::GhostClustersT ghostClusters,
    TaskHaloEngine& engine,
    const seissol::parallel::Pinning* pinning,
    bool thread)
    : AbstractCommunicationManager(std::move(ghostClusters)), engine_(engine) {
#ifdef SEISSOL_HOST_TASK_GRAPH
  if (thread) {
    helper_ = std::make_unique<seissol::parallel::HelperThread>(
        [this]() {
          engine_.poll();
          // done once the main thread has created all exchanges of the epoch, and all have
          // completed
          return allCreated_.load(std::memory_order_acquire) && engine_.outstanding() == 0;
        },
        pinning);
  }
#endif
}

void TaskCommunicationManager::progression() {
  // the ghost clusters create their tasks on the main thread
  poll();
#ifdef SEISSOL_HOST_TASK_GRAPH
  if (helper_ == nullptr) {
    engine_.poll();
  }
#endif
}

bool TaskCommunicationManager::checkIfFinished() const {
  return std::all_of(ghostClusters_.begin(), ghostClusters_.end(), [](const auto& ghostCluster) {
    return ghostCluster->synced();
  });
}

void TaskCommunicationManager::reset(double newSyncTime) {
  if (helper_ != nullptr) {
    helper_->stop();
  }
  allCreated_.store(false, std::memory_order_release);
  AbstractCommunicationManager::reset(newSyncTime);
  if (helper_ != nullptr) {
    helper_->start();
  }
}

void TaskCommunicationManager::awaitExchanges() {
  allCreated_.store(true, std::memory_order_release);
#ifdef SEISSOL_HOST_TASK_GRAPH
  if (helper_ == nullptr) {
    // the exchanges of the other processes may wait for those of this one
    while (engine_.outstanding() > 0) {
      engine_.poll();
#pragma omp taskyield
    }
  }
#endif
}

} // namespace seissol::solver

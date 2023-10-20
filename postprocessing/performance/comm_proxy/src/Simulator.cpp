// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "Simulator.hpp"

#include <Kernels/precision.hpp>
#include <Parallel/MPI.h>
#include <Solver/time_stepping/DirectGhostTimeCluster.h>

#include <algorithm>
#include <cstdint>
#include <utility>

namespace seissol::time_stepping {
class DummyTimeCluster : public AbstractTimeCluster {
  public:
  DummyTimeCluster(double maxTimeStepSize, long timeStepRate)
      : AbstractTimeCluster(maxTimeStepSize, timeStepRate) {}

  protected:
  void start() override {}
  void predict() override {}
  void correct() override {}
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override {}
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override {}
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override {}
};
} // namespace seissol::time_stepping

namespace seissol {

Simulator::Simulator(std::vector<TimeCluster> const& ltsStructure, std::shared_ptr<Allocator> alloc)
    : alloc_(std::move(alloc)) {
  auto const initMeshStructure = [&]() {
    meshStructure_.reserve(ltsStructure.size());

    std::size_t bufferSize = 0;
    for (auto const& tc : ltsStructure) {
      for (auto const& reg : tc.regions) {
        bufferSize += reg.ghostRegionSize;
        bufferSize += reg.copyRegionSize;
      }
      bufferSize += tc.interiorLayerSize;
    }
    buffer_ = (real*)alloc_->malloc(bufferSize * sizeof(real));

    real* bufferPtr = buffer_;
    for (auto const& tc : ltsStructure) {
      meshStructure_.push_back(MeshStructure{});
      auto& mc = meshStructure_.back();
      unsigned numberOfRegions = tc.regions.size();
      mc.neighboringClusters = new int[numberOfRegions][2];
      mc.numberOfGhostRegionCells = new unsigned[numberOfRegions];
      mc.numberOfGhostRegionDerivatives = new unsigned[numberOfRegions];
      mc.ghostRegions = new real*[numberOfRegions];
      mc.ghostRegionSizes = new unsigned[numberOfRegions];
      mc.numberOfCopyRegionCells = new unsigned[numberOfRegions];
      mc.numberOfCommunicatedCopyRegionDerivatives = new unsigned[numberOfRegions];
      mc.copyRegions = new real*[numberOfRegions];
      mc.copyRegionSizes = new unsigned[numberOfRegions];
      mc.sendIdentifiers = new int[numberOfRegions];
      mc.receiveIdentifiers = new int[numberOfRegions];
      mc.sendRequests = new MPI_Request[numberOfRegions];
      mc.receiveRequests = new MPI_Request[numberOfRegions];

      mc.numberOfRegions = numberOfRegions;
      mc.numberOfGhostCells = 0;
      mc.numberOfCopyCells = 0;
      mc.numberOfInteriorCells = tc.numberOfInteriorCells;
      for (unsigned r = 0; r < numberOfRegions; ++r) {
        auto& reg = tc.regions[r];
        mc.neighboringClusters[r][0] = reg.neighborRank;
        mc.neighboringClusters[r][1] = reg.neighborTimeClusterId;
        mc.numberOfGhostCells += reg.numberOfGhostRegionCells;
        mc.numberOfGhostRegionCells[r] = reg.numberOfGhostRegionCells;
        mc.numberOfGhostRegionDerivatives[r] = reg.numberOfGhostRegionDerivatives;
        mc.ghostRegions[r] = bufferPtr;
        bufferPtr += reg.ghostRegionSize;
        mc.ghostRegionSizes[r] = reg.ghostRegionSize;
        mc.numberOfCopyCells += reg.numberOfCopyRegionCells;
        mc.numberOfCopyRegionCells[r] = reg.numberOfCopyRegionCells;
        mc.numberOfCommunicatedCopyRegionDerivatives[r] =
            reg.numberOfCommunicatedCopyRegionDerivatives;
        mc.copyRegions[r] = bufferPtr;
        bufferPtr += reg.copyRegionSize;
        mc.copyRegionSizes[r] = reg.copyRegionSize;
        mc.sendIdentifiers[r] = reg.sendIdentifier;
        mc.receiveIdentifiers[r] = reg.receiveIdentifier;
        mc.sendRequests[r] = MPI_REQUEST_NULL;
        mc.receiveRequests[r] = MPI_REQUEST_NULL;
      }
      bufferPtr += tc.interiorLayerSize;
    }
  };

  auto const addClusters = [&]() {
    for (unsigned int localClusterId = 0; localClusterId < ltsStructure.size(); localClusterId++) {
      auto& tc = ltsStructure[localClusterId];
      auto* meshStructure = &meshStructure_[localClusterId];

      for (int i = 0; i < 2; ++i) {
        clusters_.push_back(std::make_unique<time_stepping::DummyTimeCluster>(tc.cflTimeStepWidth,
                                                                              tc.timeStepRate));
      }
      auto& interior = clusters_[clusters_.size() - 1];
      auto& copy = clusters_[clusters_.size() - 2];

      interior->connect(*copy);

      if (localClusterId > 0) {
        for (int i = 0; i < 2; ++i) {
          copy->connect(*clusters_[clusters_.size() - 2 - i - 1]);
          interior->connect(*clusters_[clusters_.size() - 2 - i - 1]);
        }
      }

      auto regions = std::vector<Region>(tc.regions.begin(), tc.regions.end());
      std::sort(regions.begin(), regions.end(), [](Region const& a, Region const& b) {
        return a.neighborTimeClusterId < b.neighborTimeClusterId;
      });
      auto last = std::unique(regions.begin(), regions.end(), [](Region const& a, Region const& b) {
        return a.neighborTimeClusterId == b.neighborTimeClusterId;
      });
      regions.erase(last, regions.end());

      for (auto const& reg : regions) {
        ghostClusters_.push_back(
            std::make_unique<time_stepping::DirectGhostTimeCluster>(reg.neighborCflTimeStepWidth,
                                                                    reg.neighborTimeStepRate,
                                                                    tc.timeClusterId,
                                                                    reg.neighborTimeClusterId,
                                                                    meshStructure));

        ghostClusters_.back()->connect(*copy);
      }
    }
  };

  initMeshStructure();
  addClusters();
}

Simulator::~Simulator() {
  for (auto const& mc : meshStructure_) {
    delete[] mc.neighboringClusters;
    delete[] mc.numberOfGhostRegionCells;
    delete[] mc.numberOfGhostRegionDerivatives;
    delete[] mc.ghostRegions;
    delete[] mc.ghostRegionSizes;
    delete[] mc.numberOfCopyRegionCells;
    delete[] mc.numberOfCommunicatedCopyRegionDerivatives;
    delete[] mc.copyRegions;
    delete[] mc.copyRegionSizes;
    delete[] mc.sendIdentifiers;
    delete[] mc.receiveIdentifiers;
    delete[] mc.sendRequests;
    delete[] mc.receiveRequests;
  }
  alloc_->free(buffer_);
}

void Simulator::simulate(double synchronizationTime) {
  for (auto& cluster : clusters_) {
    cluster->setSyncTime(synchronizationTime);
    cluster->reset();
  }

  for (auto& ghostCluster : ghostClusters_) {
    ghostCluster->setSyncTime(synchronizationTime);
    ghostCluster->reset();
  }

  seissol::MPI::mpi.barrier(seissol::MPI::mpi.comm());

  for (auto& cluster : clusters_) {
    cluster->act();
  }

  int maxsteps = 10;
  bool finished = false;
  while (!finished /* && maxsteps--*/) {
    finished = true;

    for (auto& cluster : clusters_) {
      poll();
      cluster->act();
    }

    finished = std::all_of(clusters_.begin(), clusters_.end(), [](auto& c) { return c->synced(); });
    finished &= poll();
  }
}

bool Simulator::poll() {
  bool finished = true;
  for (auto& ghostCluster : ghostClusters_) {
    ghostCluster->act();
    finished = finished && ghostCluster->synced();
  }
  return finished;
}

} // namespace seissol

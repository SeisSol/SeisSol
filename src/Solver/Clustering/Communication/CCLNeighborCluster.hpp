#pragma once

#include "NeighborCluster.hpp"
#include <Parallel/Runtime/Stream.hpp>
#include <atomic>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class CCLSendNeighborCluster : public SendNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  CCLSendNeighborCluster(const std::vector<RemoteCluster>& remote, void* comm);
  ~CCLSendNeighborCluster() override;

  private:
  std::vector<RemoteCluster> remote;
  std::vector<void*> memoryHandles;
  void* comm;
};

class CCLRecvNeighborCluster : public RecvNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  CCLRecvNeighborCluster(const std::vector<RemoteCluster>& remote, void* comm);
  ~CCLRecvNeighborCluster() override;

  private:
  std::vector<RemoteCluster> remote;
  std::vector<void*> memoryHandles;
  void* comm;
};

} // namespace seissol::solver::clustering::communication

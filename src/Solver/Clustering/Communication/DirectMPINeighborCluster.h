#pragma once

#include "NeighborCluster.h"
#include <Parallel/Runtime/Stream.h>
#include <atomic>
#include <mpi.h>
#include <vector>

namespace seissol::solver::clustering::communication {

class DirectMPISendNeighborCluster : public SendNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  DirectMPISendNeighborCluster(const std::vector<RemoteCluster>& remote);
  ~DirectMPISendNeighborCluster() override;

  private:
  std::vector<MPI_Request> requests;
  std::vector<int> status;
  std::mutex requestMutex;
};

class DirectMPIRecvNeighborCluster : public RecvNeighborCluster {
  public:
  bool poll() override;
  void start(parallel::runtime::StreamRuntime& runtime) override;
  void stop(parallel::runtime::StreamRuntime& runtime) override;

  DirectMPIRecvNeighborCluster(const std::vector<RemoteCluster>& remote);
  ~DirectMPIRecvNeighborCluster() override;

  private:
  std::vector<MPI_Request> requests;
  std::vector<int> status;
  std::mutex requestMutex;
};

/*
class CopyMPISendNeighborCluster : public SendNeighborCluster {
public:
    CopyMPISendNeighborCluster(const std::vector<RemoteCluster>& remote);
    ~CopyMPISendNeighborCluster() override;

    bool poll() override;
    void start(parallel::runtime::StreamRuntime& runtime) override;
private:
    std::vector<MPI_Request> requests;
    std::vector<void*> copyStreams;
    std::mutex requestMutex;
};

class CopyMPIRecvNeighborCluster : public RecvNeighborCluster {
public:
    CopyMPIRecvNeighborCluster(const std::vector<RemoteCluster>& remote);
    ~CopyMPIRecvNeighborCluster() override;

    bool poll() override;
    void start(parallel::runtime::StreamRuntime& runtime) override;
private:
    std::vector<MPI_Request> requests;
    std::vector<void*> copyStreams;
    std::mutex requestMutex;
};
*/

} // namespace seissol::solver::clustering::communication

#pragma once

#include <Initializer/Typedefs.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/Clustering/ActorState.h>
namespace seissol::solver::clustering::communication {

struct RemoteCluster {
  void* data;
  std::size_t size;
  MPI_Datatype datatype;
  int rank;
  int tag;
};

struct HaloCommunication {
  std::vector<std::vector<RemoteCluster>> ghost;
  std::vector<std::vector<RemoteCluster>> copy;
};

HaloCommunication getHaloCommunication(std::size_t clusterCount, const MeshStructure* structure);

struct CommunicationSetup {
  ComputeStep start;
  ComputeStep step;
};

class NeighborCluster {
  public:
  virtual bool poll() = 0;
  virtual void start(parallel::runtime::StreamRuntime& runtime) = 0;
  virtual void stop(parallel::runtime::StreamRuntime& runtime) = 0;
  virtual ~NeighborCluster() = default;

  void startFrom(parallel::runtime::StreamRuntime& runtime);
  void stopTo(parallel::runtime::StreamRuntime& runtime);

  private:
  parallel::runtime::StreamRuntime myRuntime;
};

class SendNeighborCluster : public NeighborCluster {};

class RecvNeighborCluster : public NeighborCluster {};

} // namespace seissol::solver::clustering::communication

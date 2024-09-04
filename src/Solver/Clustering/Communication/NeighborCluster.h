#pragma once

#include <Initializer/typedefs.hpp>
#include <Parallel/Runtime/Stream.hpp>
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

class NeighborCluster {
  public:
  virtual bool poll() = 0;
  virtual void start(parallel::runtime::StreamRuntime& runtime) = 0;
  virtual void stop(parallel::runtime::StreamRuntime& runtime) = 0;
  virtual ~NeighborCluster() = default;
};

class SendNeighborCluster : public NeighborCluster {};

class RecvNeighborCluster : public NeighborCluster {};

} // namespace seissol::solver::clustering::communication

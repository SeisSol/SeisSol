// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTER_H_

#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Initializer/Typedefs.h>
#include <Parallel/Host/CpuExecutor.h>
#include <Parallel/Runtime/Stream.h>
#include <Solver/Clustering/ActorState.h>

#include "Parallel/MPI.h"

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

HaloCommunication getHaloCommunication(const initializer::ClusterLayout& layout,
                                       const MeshStructure* structure);

struct CommunicationSetup {
  ComputeStep start;
  ComputeStep step;
};

class NeighborCluster {
  public:
  NeighborCluster(const std::shared_ptr<parallel::host::CpuExecutor>&, double priority);
  virtual bool poll() = 0;
  virtual void start(parallel::runtime::StreamRuntime& runtime) = 0;
  virtual void stop(parallel::runtime::StreamRuntime& runtime) = 0;
  virtual ~NeighborCluster() = default;

  void startFrom(parallel::runtime::StreamRuntime& runtime);
  void stopTo(parallel::runtime::StreamRuntime& runtime);
  bool checkDone();

  private:
  parallel::runtime::StreamRuntime myRuntime;
};

class SendNeighborCluster : public NeighborCluster {
  public:
  SendNeighborCluster(const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                      double priority)
      : NeighborCluster(cpuExecutor, priority) {}
  ~SendNeighborCluster() override = default;
};

class RecvNeighborCluster : public NeighborCluster {
  public:
  RecvNeighborCluster(const std::shared_ptr<parallel::host::CpuExecutor>& cpuExecutor,
                      double priority)
      : NeighborCluster(cpuExecutor, priority) {}
  ~RecvNeighborCluster() override = default;
};

} // namespace seissol::solver::clustering::communication
#endif // SEISSOL_SRC_SOLVER_CLUSTERING_COMMUNICATION_NEIGHBORCLUSTER_H_

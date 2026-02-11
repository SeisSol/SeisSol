// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_SHMEMCLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_SHMEMCLUSTER_H_

#include "Memory/MemoryAllocator.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

namespace seissol::time_stepping {

class ShmemCluster : public AbstractTimeCluster {
  protected:
  bool mayCorrect() override;
  void start() override;
  void predict() override;
  void correct() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  private:
  std::size_t globalClusterId;
  std::size_t otherGlobalClusterId;

  solver::RemoteClusterPair meshStructure;
  std::shared_ptr<parallel::runtime::StreamRuntime> stream;

  std::vector<void*> remotePointers;
  memory::MemkindArray<uint64_t> signal;

#ifdef ACL_DEVICE
  device::DeviceGraphHandle handle;
#endif

  void launchSend();

  public:
  ShmemCluster(double maxTimeStepSize,
               int timeStepRate,
               int globalTimeClusterId,
               int otherGlobalTimeClusterId,
               const seissol::solver::HaloCommunication& meshStructure,
               bool persistent);
  ~ShmemCluster() override;

  [[nodiscard]] std::string description() const override;
};

} // namespace seissol::time_stepping
#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_SHMEMCLUSTER_H_

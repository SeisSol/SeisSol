// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLNEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLNEIGHBORCLUSTER_H_

#include <Parallel/Runtime/Stream.h>
#include <Solver/TimeStepping/AbstractTimeCluster.h>
#include <Solver/TimeStepping/HaloCommunication.h>
namespace seissol::time_stepping {

class CCLNeighborCluster : public AbstractTimeCluster {
  protected:
  bool mayCorrect() override;
  void start() override;
  void predict() override;
  void correct() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  private:
  std::shared_ptr<parallel::runtime::StreamRuntime> stream;
  std::vector<solver::RemoteCluster> remote;
  std::vector<bool> isSend;
  std::vector<void*> memoryHandles;
  void* comm;
  void* event{nullptr};

#ifdef ACL_DEVICE
  device::DeviceGraphHandle handle;
#endif

  public:
  CCLNeighborCluster(double maxTimeStepSize,
                     int timeStepRate,
                     int globalTimeClusterId,
                     int otherGlobalTimeClusterId,
                     const seissol::solver::HaloCommunication& meshStructure,
                     bool persistent);
  ~CCLNeighborCluster() override;
};

} // namespace seissol::time_stepping
#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLNEIGHBORCLUSTER_H_

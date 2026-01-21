// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLNEIGHBORCLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLNEIGHBORCLUSTER_H_

#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

namespace seissol::time_stepping {

std::vector<void*> createComms(std::size_t count);

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
  void* commSend;
  void* commRecv;
  void* event{nullptr};

#ifdef ACL_DEVICE
  device::DeviceGraphHandle handle1;
  device::DeviceGraphHandle handle2;
  device::DeviceGraphHandle handle3;
#endif

  void launch(bool send, bool recv);

  std::size_t globalClusterId;
  std::size_t otherGlobalClusterId;

  public:
  CCLNeighborCluster(double maxTimeStepSize,
                     int timeStepRate,
                     int globalTimeClusterId,
                     int otherGlobalTimeClusterId,
                     const seissol::solver::HaloCommunication& meshStructure,
                     bool persistent,
                     const std::vector<void*>& comms);
  ~CCLNeighborCluster() override;

  [[nodiscard]] std::string description() const override;
};

} // namespace seissol::time_stepping
#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_CCLNEIGHBORCLUSTER_H_

// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTGHOSTTIMECLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTGHOSTTIMECLUSTER_H_

#include "AbstractTimeCluster.h"
#include "Initializer/Typedefs.h"
#include <list>
#include <string>

#include "Parallel/MPI.h"

#include "Solver/TimeStepping/HaloCommunication.h"

namespace seissol::time_stepping {
class AbstractGhostTimeCluster : public AbstractTimeCluster {
  protected:
  std::size_t globalClusterId;
  std::size_t otherGlobalClusterId;
  solver::RemoteClusterPair meshStructure;
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;
  std::list<unsigned int> sendQueue;
  std::list<unsigned int> receiveQueue;

  double lastSendTime = -1.0;

  virtual void sendCopyLayer() = 0;
  virtual void receiveGhostLayer() = 0;

  bool testQueue(MPI_Request* requests, std::list<unsigned int>& regions);
  bool testForCopyLayerSends();
  virtual bool testForGhostLayerReceives() = 0;

  void start() override;
  void predict() override;
  void correct() override;
  bool mayPredict() override;
  bool mayCorrect() override;
  bool maySync() override;
  void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
  void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  public:
  AbstractGhostTimeCluster(double maxTimeStepSize,
                           int timeStepRate,
                           int globalTimeClusterId,
                           int otherGlobalTimeClusterId,
                           const seissol::solver::HaloCommunication& meshStructure);

  void reset() override;
  ActResult act() override;

  [[nodiscard]] std::string description() const override;
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTGHOSTTIMECLUSTER_H_

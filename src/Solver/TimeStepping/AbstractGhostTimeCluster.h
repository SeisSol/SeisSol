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
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <list>
#include <string>

namespace seissol::time_stepping {
class AbstractGhostTimeCluster : public AbstractTimeCluster {
  protected:
  std::size_t globalClusterId;
  std::size_t otherGlobalClusterId;
  solver::RemoteClusterPair meshStructure;
  std::vector<MPI_Request> sendRequests;
  std::vector<MPI_Request> recvRequests;
  std::list<std::size_t> sendQueue;
  std::list<std::size_t> receiveQueue;

  std::string displayName_;
  std::string otherDisplayName_;

  double lastSendTime = -1.0;

  virtual void sendCopyLayer() = 0;
  virtual void receiveGhostLayer() = 0;

  bool testQueue(MPI_Request* requests, std::list<std::size_t>& regions);
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

  [[nodiscard]] bool timeoutFail() const override;

  public:
  AbstractGhostTimeCluster(double maxTimeStepSize,
                           std::uint64_t timeStepRate,
                           std::size_t globalTimeClusterId,
                           std::size_t otherGlobalTimeClusterId,
                           const std::string& displayName,
                           const std::string& otherDisplayName,
                           const seissol::solver::HaloCommunication& meshStructure);

  void reset() override;
  ActResult act() override;

  [[nodiscard]] std::string description() const override;
};
} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_ABSTRACTGHOSTTIMECLUSTER_H_

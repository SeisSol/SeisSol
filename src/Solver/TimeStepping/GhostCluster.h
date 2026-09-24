// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTCLUSTER_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTCLUSTER_H_

#include "AbstractTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"
#include "Solver/TimeStepping/HaloTransport.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace seissol::time_stepping {

/**
 * Stands in for one remote time cluster next to the copy layer of a local one: it follows the
 * progress of the copy layer, decides when the halo data is exchanged, and lets the transport move
 * it.
 */
class GhostCluster : public AbstractTimeCluster {
  public:
  GhostCluster(double maxTimeStepSize,
               std::uint64_t timeStepRate,
               const std::string& displayName,
               const std::string& otherDisplayName,
               const solver::RemoteClusterPair& regions,
               std::unique_ptr<HaloTransport> transport);

  void reset() override;
  ActResult act() override;
  void finalize() override;

  [[nodiscard]] std::string description() const override;

  /// Number of messages sent so far, one per copy region and exchange.
  [[nodiscard]] std::size_t sentMessages() const;

  /// Number of messages received so far, one per ghost region and exchange.
  [[nodiscard]] std::size_t receivedMessages() const;

  protected:
  void start() override;
  void predict() override;
  void correct() override;
  bool mayPredict() override;
  bool mayCorrect() override;
  bool maySync() override;
  void handleNeighborPrediction(const NeighborCluster& neighbor) override;
  void handleNeighborCorrection(const NeighborCluster& neighbor) override;

  [[nodiscard]] bool timeoutFail() const override;

  private:
  void sendCopyLayer();
  void receiveGhostLayer();
  bool testForCopyLayerSends();
  bool testForGhostLayerReceives();

  std::unique_ptr<HaloTransport> transport_;
  std::size_t copyRegionCount_;
  std::size_t ghostRegionCount_;
  std::size_t sentMessages_{0};
  std::size_t receivedMessages_{0};

  std::string displayName_;
  std::string otherDisplayName_;

  double lastSendTime_ = -1.0;
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTCLUSTER_H_

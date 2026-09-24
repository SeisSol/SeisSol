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

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

namespace seissol::time_stepping {

/**
 * Stands in for one remote time cluster next to the copy layer of a local one: it follows the
 * progress of the copy layer, decides when the halo data is exchanged, and lets the transport move
 * it.
 *
 * Receiving and sending progress independently of each other. The progress of the receives is
 * published as the predictions of this cluster: the ghost data is in place up to there, and the
 * copy layer may correct up to there. The progress of the sends is published as its corrections:
 * the copy data is out up to there, and the copy layer may predict beyond. Both advance in steps
 * of the exchange period, the larger one of the two time step rates, and end with the last step
 * of this cluster before the synchronization point.
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
  void predict() override {}
  void correct() override {}
  void handleNeighborPrediction(const NeighborCluster& neighbor) override;
  void handleNeighborCorrection(const NeighborCluster& neighbor) override;

  [[nodiscard]] bool timeoutFail() const override;
  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

  private:
  void sendCopyLayer(long target);
  void receiveGhostLayer(long target);
  bool testForCopyLayerSends();
  bool testForGhostLayerReceives();

  /// Moves the receive progress (and the prediction time) forward to `target`.
  void advanceReceived(long target);

  /// Moves the send progress (and the correction time) forward to `target`.
  void advanceSent(long target);

  /// The progress at the synchronization point: the end of the last step before it.
  [[nodiscard]] long finalSteps() const;

  /// The number of steps between two exchanges.
  [[nodiscard]] long exchangePeriod() const;

  std::unique_ptr<HaloTransport> transport_;
  std::size_t copyRegionCount_;
  std::size_t ghostRegionCount_;
  std::size_t sentMessages_{0};
  std::size_t receivedMessages_{0};

  std::string displayName_;
  std::string otherDisplayName_;

  bool receiving_{false};
  bool sending_{false};
  long receiveTarget_{0};
  long sendTarget_{0};
};

} // namespace seissol::time_stepping

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_GHOSTCLUSTER_H_

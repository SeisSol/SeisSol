// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_SPLITHALOTRANSPORT_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_SPLITHALOTRANSPORT_H_

#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::solver {

/**
 * Moves the regions with the processes of the own node with one transport (e.g. CCL, ordered on the
 * device), and those with the processes of other nodes with another one, which the host drives
 * (e.g. GPU-aware MPI): a send or receive of the latter starts once the work behind the event it
 * was started after has completed on the device, as the host sees by polling. The transport is
 * ordered on the device as the first one is; a send or receive counts as done once the second one
 * has completed it, and the first one has launched (or completed) it.
 *
 * Experimental: so that the latency of the network does not add up along the order of the
 * exchanges on the device, while the exchange within a node stays on the device.
 */
class SplitHaloTransport : public HaloTransport {
  public:
  SplitHaloTransport(std::unique_ptr<HaloTransport> local, std::unique_ptr<HaloTransport> remote);
  ~SplitHaloTransport() override = default;

  SplitHaloTransport(const SplitHaloTransport&) = delete;
  SplitHaloTransport(SplitHaloTransport&&) = delete;
  SplitHaloTransport& operator=(const SplitHaloTransport&) = delete;
  SplitHaloTransport& operator=(SplitHaloTransport&&) = delete;

  /**
   * Splits the regions into those with the processes `local` says are on the own node, and the
   * others.
   */
  static std::pair<RemoteClusterPair, RemoteClusterPair> split(const RemoteClusterPair& regions,
                                                               const std::vector<bool>& local);

  void startSend() override;
  bool testSend() override;
  void startReceive() override;
  bool testReceive() override;

  [[nodiscard]] bool streamOrdered() const override { return local_->streamOrdered(); }
  [[nodiscard]] std::string describe() const override;
  [[nodiscard]] ActorEvent latestEvent() const override { return local_->latestEvent(); }

  void startSendAfter(const ActorEvent& event) override;
  void startReceiveAfter(const ActorEvent& event) override;
  void startInterval(const ExchangeInterval& interval) override;
  void finalize() override;
  void setTraceName(const std::string& name) override;

  private:
  /// starts the deferred operations of the remote transport whose events have completed
  void startDeferred();

  std::unique_ptr<HaloTransport> local_;
  std::unique_ptr<HaloTransport> remote_;
  bool sendDeferred_{false};
  bool receiveDeferred_{false};
  ActorEvent sendAfter_;
  ActorEvent receiveAfter_;
};

} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALO_SPLITHALOTRANSPORT_H_

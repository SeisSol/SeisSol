// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Executor.h"
#include "Kernels/Precision.h"
#include "Solver/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Solver/TimeStepping/Actor/ActorState.h"
#include "Solver/TimeStepping/Halo/GhostCluster.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::solver;

namespace {

/**
 * A transport that completes each send and receive after a given number of tests.
 */
class FakeTransport : public HaloTransport {
  public:
  FakeTransport(int sendLatency, int receiveLatency)
      : sendLatency_(sendLatency), receiveLatency_(receiveLatency) {}

  void startSend() override {
    REQUIRE_FALSE(sending);
    sending = true;
    sendTests_ = 0;
  }
  bool testSend() override {
    if (sending && ++sendTests_ > sendLatency_) {
      sending = false;
      ++completedSends;
    }
    return !sending;
  }
  void startReceive() override {
    REQUIRE_FALSE(receiving);
    receiving = true;
    receiveTests_ = 0;
  }
  bool testReceive() override {
    if (receiving && ++receiveTests_ > receiveLatency_) {
      receiving = false;
      ++completedReceives;
    }
    return !receiving;
  }

  bool sending{false};
  bool receiving{false};
  long completedSends{0};
  long completedReceives{0};

  private:
  int sendLatency_;
  int receiveLatency_;
  int sendTests_{0};
  int receiveTests_{0};
};

/**
 * A copy layer which checks, whenever it acts, that the halo exchange has delivered what the step
 * needs: the ghost data up to the end of a corrected step, and the copy data up to the start of a
 * predicted step.
 */
class CopyCluster : public AbstractTimeCluster {
  public:
  CopyCluster(long timeStepRate, long exchangePeriod, const FakeTransport& transport)
      : AbstractTimeCluster(static_cast<double>(timeStepRate), timeStepRate, Executor::Host),
        exchangePeriod_(exchangePeriod), transport_(transport) {}

  void setInterval(long receivesBefore, long sendsBefore) {
    receivesBefore_ = receivesBefore;
    sendsBefore_ = sendsBefore;
  }

  long correctionsDuringSends{0};

  protected:
  void start() override {}
  void predict() override {
    // the copy data of all complete exchanges before this step has been sent
    const auto start = ct_.predictionsSinceLastSync;
    CHECK(transport_.completedSends - sendsBefore_ >= start / exchangePeriod_);
  }
  void correct() override {
    // the ghost data up to the end of this step has arrived
    const auto end = std::min(ct_.stepsSinceLastSync + ct_.timeStepRate, ct_.stepsUntilSync);
    CHECK(transport_.completedReceives - receivesBefore_ >=
          (end + exchangePeriod_ - 1) / exchangePeriod_);
    if (transport_.sending) {
      ++correctionsDuringSends;
    }
  }
  void handleNeighborPrediction(const NeighborCluster& /*neighbor*/) override {}
  void handleNeighborCorrection(const NeighborCluster& /*neighbor*/) override {}
  void printTimeoutMessage(std::chrono::seconds /*timeSinceLastUpdate*/) override {}

  private:
  long exchangePeriod_;
  const FakeTransport& transport_;
  long receivesBefore_{0};
  long sendsBefore_{0};
};

solver::RemoteClusterPair oneRegionEach() {
  solver::RemoteClusterPair regions;
  regions.copy.emplace_back(nullptr, 1, RealType::F64, 1, 0);
  regions.ghost.emplace_back(nullptr, 1, RealType::F64, 1, 0);
  return regions;
}

long roundUp(long value, long multiple) { return (value + multiple - 1) / multiple * multiple; }

long divideUp(long value, long divisor) { return (value + divisor - 1) / divisor; }

/**
 * Runs a copy layer next to a ghost cluster through the given synchronization points, and checks
 * after each interval that exactly as many sends and receives have happened as the remote side
 * expects. Returns the number of copy corrections that happened while a send was in flight.
 */
long exchange(long copyRate,
              long ghostRate,
              int sendLatency,
              int receiveLatency,
              bool ghostFirst,
              const std::vector<long>& syncTimes) {
  const auto exchangePeriod = std::max(copyRate, ghostRate);
  auto transportOwner = std::make_unique<FakeTransport>(sendLatency, receiveLatency);
  auto& transport = *transportOwner;
  CopyCluster copy(copyRate, exchangePeriod, transport);
  GhostCluster ghost(static_cast<double>(ghostRate),
                     ghostRate,
                     "copy",
                     "ghost",
                     oneRegionEach(),
                     std::move(transportOwner));
  copy.connect(ghost);

  long start = 0;
  for (const auto syncTime : syncTimes) {
    copy.setInterval(transport.completedReceives, transport.completedSends);
    const auto sendsBefore = transport.completedSends;
    const auto receivesBefore = transport.completedReceives;

    copy.setSyncTime(static_cast<double>(syncTime));
    ghost.setSyncTime(static_cast<double>(syncTime));
    copy.reset();
    ghost.reset();

    std::size_t iterations = 0;
    while (!copy.synced() || !ghost.synced() || iterations == 0) {
      REQUIRE(iterations < 100000);
      ++iterations;
      if (ghostFirst) {
        ghost.act();
        copy.act();
      } else {
        copy.act();
        ghost.act();
      }
    }

    // the remote side sends once per exchange up to the end of the last step of the ghost cluster,
    // and receives once per exchange up to the end of the last step of the copy layer
    const auto steps = syncTime - start;
    CHECK(transport.completedReceives - receivesBefore ==
          divideUp(roundUp(steps, ghostRate), exchangePeriod));
    CHECK(transport.completedSends - sendsBefore ==
          divideUp(roundUp(steps, copyRate), exchangePeriod));
    CHECK_FALSE(transport.sending);
    CHECK_FALSE(transport.receiving);
    start = syncTime;
  }
  CHECK(ghost.sentMessages() == static_cast<std::size_t>(transport.completedSends));
  CHECK(ghost.receivedMessages() == static_cast<std::size_t>(transport.completedReceives));
  return copy.correctionsDuringSends;
}

} // namespace

TEST_CASE("Ghost clusters exchange the halo data" * doctest::test_suite("solver")) {
  // synchronization points on and between the steps of the larger cluster
  const std::vector<long> syncTimes{4, 7, 12, 13, 20};

  for (const auto ghostFirst : {false, true}) {
    for (const auto [copyRate, ghostRate] :
         std::vector<std::pair<long, long>>{{1, 1}, {2, 2}, {1, 2}, {2, 1}, {1, 3}, {3, 1}}) {
      for (const auto [sendLatency, receiveLatency] :
           std::vector<std::pair<int, int>>{{0, 0}, {3, 0}, {0, 3}, {2, 5}}) {
        CAPTURE(ghostFirst);
        CAPTURE(copyRate);
        CAPTURE(ghostRate);
        CAPTURE(sendLatency);
        CAPTURE(receiveLatency);
        exchange(copyRate, ghostRate, sendLatency, receiveLatency, ghostFirst, syncTimes);
      }
    }
  }
}

TEST_CASE("Ghost data arrives while the copy data is still being sent" *
          doctest::test_suite("solver")) {
  // Next to a smaller remote cluster, the copy layer receives and sends once per own step. A slow
  // send of that step must not hold back its correction.
  CHECK(exchange(2, 1, 5, 0, false, {20}) > 0);
  CHECK(exchange(3, 1, 5, 0, true, {20}) > 0);
}

} // namespace seissol::unit_test

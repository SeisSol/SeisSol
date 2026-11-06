// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Common/Real.h"
#include "Solver/TimeStepping/Halo/HaloCheck.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::solver;

namespace halocheck {

constexpr int Rank = 3;
constexpr std::size_t TagF64 = 11;
constexpr std::size_t TagF32 = 12;

/**
 * The regions of a process that exchanges with itself: the copy region of each tag goes to the
 * ghost region of the same tag. They start out with values that the check never writes.
 */
struct Loopback {
  std::vector<double> copyF64;
  std::vector<double> ghostF64;
  std::vector<float> copyF32;
  std::vector<float> ghostF32;
  RemoteClusterPair regions;

  /// `sender` is where the ghost regions expect their data from
  Loopback(std::size_t sizeF64, std::size_t sizeF32, int sender = Rank)
      : copyF64(sizeF64, -1.0), ghostF64(sizeF64, -2.0), copyF32(sizeF32, -3.0F),
        ghostF32(sizeF32, -4.0F) {
    regions.copy.emplace_back(copyF64.data(), sizeF64, RealType::F64, Rank, TagF64);
    regions.copy.emplace_back(copyF32.data(), sizeF32, RealType::F32, Rank, TagF32);
    regions.ghost.emplace_back(ghostF64.data(), sizeF64, RealType::F64, sender, TagF64);
    regions.ghost.emplace_back(ghostF32.data(), sizeF32, RealType::F32, sender, TagF32);
  }

  /// moves the copy regions into the ghost regions, as the network would; without the last
  /// `missing` values of the first one
  void deliver(std::size_t missing) {
    std::copy(
        copyF64.begin(), copyF64.end() - static_cast<std::ptrdiff_t>(missing), ghostF64.begin());
    std::ranges::copy(copyF32, ghostF32.begin());
  }

  [[nodiscard]] bool allZero() const {
    const auto zero = [](auto value) { return value == 0; };
    return std::ranges::all_of(copyF64, zero) && std::ranges::all_of(ghostF64, zero) &&
           std::ranges::all_of(copyF32, zero) && std::ranges::all_of(ghostF32, zero);
  }
};

/**
 * Moves the data of a loopback once both its send and its receive have started; and records the
 * intervals announced to it.
 */
class LoopbackTransport : public HaloTransport {
  public:
  LoopbackTransport(Loopback& loopback, std::size_t missing)
      : loopback_(loopback), missing_(missing) {}

  void startSend() override { sending_ = true; }
  bool testSend() override {
    progress();
    return !sending_;
  }
  void startReceive() override { receiving_ = true; }
  bool testReceive() override {
    progress();
    return !receiving_;
  }
  void startInterval(const ExchangeInterval& interval) override { intervals.push_back(interval); }

  std::vector<ExchangeInterval> intervals;

  private:
  void progress() {
    if (sending_ && receiving_) {
      loopback_.deliver(missing_);
      sending_ = false;
      receiving_ = false;
    }
  }

  Loopback& loopback_;
  std::size_t missing_;
  bool sending_{false};
  bool receiving_{false};
};

inline HaloCheckResult check(HaloTransport& transport, const Loopback& loopback) {
  return checkHaloExchange({{.transport = &transport,
                             .regions = &loopback.regions,
                             .description = "copy-0 and ghost-1"}},
                           Rank,
                           HaloMemory::Host,
                           std::chrono::seconds(0),
                           false);
}

} // namespace halocheck

TEST_CASE("The values of the halo check tell the sender, the tag and the position, and single "
          "precision holds them" *
          doctest::test_suite("solver")) {
  for (const int sender : {0, 1, 7, 4095, 1 << 20}) {
    for (const std::size_t tag : {0UL, 1UL, 12345UL}) {
      for (const std::size_t position : {0UL, 1UL, 1000UL, (1UL << 24U) + 5}) {
        const auto value = haloCheckValue(sender, tag, position);
        CHECK(value >= 1);
        CHECK(value < static_cast<double>(1UL << 23U));
        CHECK(static_cast<double>(static_cast<float>(value)) == value);
      }
    }
  }
  CHECK(haloCheckValue(2, 5, 10) != haloCheckValue(2, 5, 11));
  CHECK(haloCheckValue(2, 5, 10) != haloCheckValue(3, 5, 10));
  CHECK(haloCheckValue(2, 5, 10) != haloCheckValue(2, 6, 10));
}

TEST_CASE("The halo check passes an exchange that delivers the data, and leaves the regions zero" *
          doctest::test_suite("solver")) {
  // more values than go through the buffer for the device at once
  const std::size_t sizeF64 = (1UL << 20U) + 3;
  const std::size_t sizeF32 = 1000;
  halocheck::Loopback loopback(sizeF64, sizeF32);
  halocheck::LoopbackTransport transport(loopback, 0);
  const auto result = halocheck::check(transport, loopback);
  CHECK(result.regions == 2);
  CHECK(result.bytes == sizeF64 * sizeof(double) + sizeF32 * sizeof(float));
  CHECK(result.wrongRegions == 0);
  CHECK(result.firstWrong.empty());
  // one exchange in each direction
  REQUIRE(transport.intervals.size() == 1);
  CHECK(transport.intervals[0].sendSteps == 1);
  CHECK(transport.intervals[0].receiveSteps == 1);
  CHECK(loopback.allZero());
}

TEST_CASE("The halo check finds the values that a ghost region misses" *
          doctest::test_suite("solver")) {
  halocheck::Loopback loopback(100, 10);
  halocheck::LoopbackTransport transport(loopback, 1);
  const auto result = halocheck::check(transport, loopback);
  CHECK(result.regions == 2);
  CHECK(result.wrongRegions == 1);
  // the last value of the region with tag 11 stayed zero
  CHECK(result.firstWrong.find("between copy-0 and ghost-1") != std::string::npos);
  CHECK(result.firstWrong.find("from rank 3 with tag 11 has 1 of 100 values (f64) wrong, the first "
                               "at position 99: 0 instead of") != std::string::npos);
  CHECK(loopback.allZero());
}

TEST_CASE("The halo check finds data from another sender" * doctest::test_suite("solver")) {
  halocheck::Loopback loopback(100, 10, halocheck::Rank + 1);
  halocheck::LoopbackTransport transport(loopback, 0);
  const auto result = halocheck::check(transport, loopback);
  CHECK(result.wrongRegions == 2);
  CHECK(result.firstWrong.find("has 100 of 100 values (f64) wrong, the first at position 0") !=
        std::string::npos);
  CHECK(loopback.allZero());
}

} // namespace seissol::unit_test

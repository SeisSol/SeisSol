// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HaloCheck.h"

#include "Common/Real.h"
#include "Monitoring/Unit.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"
#include "Solver/TimeStepping/Halo/HaloTransport.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <sstream>
#include <string>
#include <utils/logger.h>
#include <vector>

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif

namespace seissol::solver {

namespace {

/// the values of a region that go through the buffer on the host at once, to and from the device
constexpr std::size_t ChunkValues = 1UL << 20U;

/**
 * Reaches the regions: directly in host memory, or through a buffer on the host that gets copied
 * to and from device memory.
 */
class RegionAccess {
  public:
  explicit RegionAccess(HaloMemory memory) : memory_(memory) {
#ifdef ACL_DEVICE
    if (memory_ == HaloMemory::Device) {
      buffer_ =
          device::DeviceInstance::instance().api().allocPinnedMem(ChunkValues * sizeof(double));
    }
#else
    assert(memory_ == HaloMemory::Host);
#endif
  }

  ~RegionAccess() {
#ifdef ACL_DEVICE
    if (buffer_ != nullptr) {
      device::DeviceInstance::instance().api().freePinnedMem(buffer_);
    }
#endif
  }

  RegionAccess(const RegionAccess&) = delete;
  RegionAccess(RegionAccess&&) = delete;
  RegionAccess& operator=(const RegionAccess&) = delete;
  RegionAccess& operator=(RegionAccess&&) = delete;

  /**
   * Writes a region: `values(position, count, target)` puts its `count` values from `position` on
   * into `target`.
   */
  template <typename T, typename F>
  void write(T* data, std::size_t size, const F& values) {
    for (std::size_t begin = 0; begin < size; begin += ChunkValues) {
      const auto count = std::min(ChunkValues, size - begin);
      if (memory_ == HaloMemory::Host) {
        values(begin, count, data + begin);
      } else {
        auto* chunk = static_cast<T*>(buffer_);
        values(begin, count, chunk);
#ifdef ACL_DEVICE
        device::DeviceInstance::instance().api().copyTo(data + begin, chunk, count * sizeof(T));
#endif
      }
    }
  }

  /**
   * Reads a region: `check(position, count, source)` gets its `count` values from `position` on in
   * `source`.
   */
  template <typename T, typename F>
  void read(const T* data, std::size_t size, const F& check) {
    for (std::size_t begin = 0; begin < size; begin += ChunkValues) {
      const auto count = std::min(ChunkValues, size - begin);
      if (memory_ == HaloMemory::Host) {
        check(begin, count, data + begin);
      } else {
        auto* chunk = static_cast<T*>(buffer_);
#ifdef ACL_DEVICE
        device::DeviceInstance::instance().api().copyFrom(chunk, data + begin, count * sizeof(T));
#endif
        check(begin, count, static_cast<const T*>(chunk));
      }
    }
  }

  private:
  HaloMemory memory_;
  void* buffer_{nullptr};
};

/// calls `function` with a value of the type of the reals of the region
template <typename F>
void withReals(const RemoteCluster& region, const F& function) {
  if (region.datatype == RealType::F32) {
    function(float{});
  } else {
    function(double{});
  }
}

void writeValues(RegionAccess& access, const RemoteCluster& region, int sender) {
  withReals(region, [&](auto real) {
    using T = decltype(real);
    access.write(static_cast<T*>(region.data),
                 region.size,
                 [&](std::size_t position, std::size_t count, T* target) {
                   for (std::size_t i = 0; i < count; ++i) {
                     target[i] = static_cast<T>(haloCheckValue(sender, region.tag, position + i));
                   }
                 });
  });
}

void writeZeros(RegionAccess& access, const RemoteCluster& region) {
  withReals(region, [&](auto real) {
    using T = decltype(real);
    access.write(static_cast<T*>(region.data),
                 region.size,
                 [](std::size_t /*position*/, std::size_t count, T* target) {
                   std::fill_n(target, count, T{0});
                 });
  });
}

/**
 * The values of a ghost region that differ from those its sender has written; what is wrong with
 * the first of them goes into `first`.
 */
std::size_t
    countWrongValues(RegionAccess& access, const RemoteCluster& region, std::string& first) {
  std::size_t wrong = 0;
  withReals(region, [&](auto real) {
    using T = decltype(real);
    access.read(static_cast<const T*>(region.data),
                region.size,
                [&](std::size_t position, std::size_t count, const T* source) {
                  for (std::size_t i = 0; i < count; ++i) {
                    const auto expected =
                        static_cast<T>(haloCheckValue(region.rank, region.tag, position + i));
                    if (source[i] != expected) {
                      if (wrong == 0) {
                        std::ostringstream text;
                        text << "at position " << position + i << ": " << source[i]
                             << " instead of " << expected;
                        first = text.str();
                      }
                      ++wrong;
                    }
                  }
                });
  });
  return wrong;
}

/**
 * Exchanges the regions of all transports once in each direction, and waits until all exchanges
 * have completed (or for the transports on device streams, have been enqueued).
 */
void exchangeOnce(const std::vector<CheckedTransport>& transports,
                  std::chrono::seconds timeout,
                  bool abortOnTimeout) {
  // the transports on device streams order the exchanges as those of an interval
  ExchangeInterval once;
  once.sendSteps = 1;
  once.receiveSteps = 1;
  for (const auto& checked : transports) {
    checked.transport->startInterval(once);
  }
  for (const auto& checked : transports) {
    checked.transport->startReceive();
  }
  for (const auto& checked : transports) {
    checked.transport->startSend();
  }

  std::vector<bool> received(transports.size(), false);
  std::vector<bool> sent(transports.size(), false);
  std::size_t pending = 2 * transports.size();
  auto deadline = std::chrono::steady_clock::now() + timeout;
  while (pending > 0) {
    for (std::size_t i = 0; i < transports.size(); ++i) {
      if (!received[i] && transports[i].transport->testReceive()) {
        received[i] = true;
        --pending;
      }
      if (!sent[i] && transports[i].transport->testSend()) {
        sent[i] = true;
        --pending;
      }
    }
    if (pending > 0 && timeout.count() > 0 && std::chrono::steady_clock::now() > deadline) {
      for (std::size_t i = 0; i < transports.size(); ++i) {
        if (!received[i] || !sent[i]) {
          std::string waits;
          if (!received[i]) {
            waits = sent[i] ? "the receive" : "the receive and the send";
          } else {
            waits = "the send";
          }
          const auto between = "between " + transports[i].description + ":";
          logWarning(true) << "The test of the halo exchange waits for" << waits.c_str()
                           << between.c_str() << transports[i].transport->describe().c_str();
        }
      }
      if (abortOnTimeout) {
        logError() << "The test of the halo exchange before the time stepping has not completed"
                   << "after" << timeout.count()
                   << "s (SEISSOL_WATCHDOG_TIMEOUT); see the warnings of the ranks above for"
                   << "what is still in flight.";
      }
      deadline += timeout;
    }
  }
}

} // namespace

double haloCheckValue(int sender, std::size_t tag, std::size_t position) {
  // the largest prime below 2^23; the factors spread the senders and the tags over it
  constexpr std::uint64_t Modulus = 8388593;
  constexpr std::uint64_t SenderFactor = 1000003;
  constexpr std::uint64_t TagFactor = 7919;
  const auto mixed = (static_cast<std::uint64_t>(sender) * SenderFactor +
                      static_cast<std::uint64_t>(tag) * TagFactor + position) %
                     Modulus;
  return static_cast<double>(mixed + 1);
}

HaloCheckResult checkHaloExchange(const std::vector<CheckedTransport>& transports,
                                  int rank,
                                  HaloMemory memory,
                                  std::chrono::seconds timeout,
                                  bool abortOnTimeout) {
  RegionAccess access(memory);

  // the ghost regions start out zero, which no value of a sender is (the copy regions do not
  // overlap: a copy cell next to several processes is there for each of them)
  for (const auto& checked : transports) {
    for (const auto& region : checked.regions->ghost) {
      writeZeros(access, region);
    }
    for (const auto& region : checked.regions->copy) {
      writeValues(access, region, rank);
    }
  }

  exchangeOnce(transports, timeout, abortOnTimeout);
#ifdef ACL_DEVICE
  if (memory == HaloMemory::Device) {
    // the transports on device streams only enqueue the exchanges, and others may receive with a
    // copy to the device
    device::DeviceInstance::instance().api().syncDevice();
  }
#endif

  HaloCheckResult result;
  for (const auto& checked : transports) {
    for (const auto& region : checked.regions->ghost) {
      ++result.regions;
      result.bytes += region.size * sizeOfRealType(region.datatype);
      std::string first;
      const auto wrong = countWrongValues(access, region, first);
      if (wrong > 0) {
        if (result.wrongRegions == 0) {
          std::ostringstream text;
          text << "between " << checked.description << ", the ghost region from rank "
               << region.rank << " with tag " << region.tag << " has " << wrong << " of "
               << region.size << " values (" << stringRealType(region.datatype)
               << ") wrong, the first " << first;
          result.firstWrong = text.str();
        }
        ++result.wrongRegions;
      }
    }
  }

  // as after the allocation
  for (const auto& checked : transports) {
    for (const auto& region : checked.regions->copy) {
      writeZeros(access, region);
    }
    for (const auto& region : checked.regions->ghost) {
      writeZeros(access, region);
    }
  }
  return result;
}

void testHaloExchange(const std::vector<CheckedTransport>& transports,
                      HaloMemory memory,
                      const std::string& mode,
                      std::chrono::seconds timeout,
                      bool abortOnTimeout) {
  std::size_t regions = 0;
  for (const auto& checked : transports) {
    regions += checked.regions->copy.size() + checked.regions->ghost.size();
  }
  MPI_Allreduce(
      MPI_IN_PLACE, &regions, 1, Mpi::castToMpiType<std::size_t>(), MPI_SUM, Mpi::mpi.comm());
  if (regions == 0) {
    return;
  }
  // one item, as the logger puts a space between its items
  const auto exchange = "The halo exchange (transfer mode " + mode + ")";
  // a crash in the library of the exchange (e.g. an MPI that cannot read device memory) comes
  // right after this
  logInfo() << "Testing the halo exchange before the time stepping, with the transfer mode"
            << mode.c_str();

  const auto start = std::chrono::steady_clock::now();
  const auto result =
      checkHaloExchange(transports, Mpi::mpi.rank(), memory, timeout, abortOnTimeout);
  double time = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  if (result.wrongRegions > 0) {
    logWarning(true) << "The halo exchange delivered wrong data to" << result.wrongRegions << "of"
                     << result.regions << "ghost regions of this rank; e.g."
                     << result.firstWrong.c_str();
  }
  std::array<std::size_t, 3> sums{result.regions, result.bytes, result.wrongRegions};
  MPI_Allreduce(MPI_IN_PLACE,
                sums.data(),
                sums.size(),
                Mpi::castToMpiType<std::size_t>(),
                MPI_SUM,
                Mpi::mpi.comm());
  MPI_Allreduce(MPI_IN_PLACE, &time, 1, MPI_DOUBLE, MPI_MAX, Mpi::mpi.comm());
  if (sums[2] > 0) {
    logError() << exchange.c_str() << "delivered wrong data to" << sums[2] << "of" << sums[0]
               << "ghost regions in the test before the time stepping; see the warnings of the"
               << "ranks above.";
  }
  logInfo() << exchange.c_str() << "works:" << sums[0] << "ghost regions received"
            << UnitByte.formatPrefix(static_cast<double>(sums[1])).c_str() << "intact (in" << time
            << "s).";
}

} // namespace seissol::solver

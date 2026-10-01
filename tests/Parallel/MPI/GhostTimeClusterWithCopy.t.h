// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// In the host transfer mode (SEISSOL_PREFERRED_MPI_DATA_TRANSFER_MODE=host), the ghost cluster
// copies the copy layer from the device to pinned host memory, sends it from there, receives the
// ghost layer into pinned host memory and copies it to the device. MPI must therefore only ever
// see the pinned buffers. Here, every rank sends two device regions around a ring, one to each
// neighbor, and checks that its device ghost regions hold exactly what the neighbors sent.

#include <doctest.h>

#ifdef ACL_DEVICE

#include "Common/Real.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/GhostTimeClusterWithCopy.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <Device/device.h>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <vector>

namespace seissol::unit_test {

namespace ghosttimeclusterwithcopytest {

// makes the steps of an exchange callable, which the actor protocol performs otherwise
class HostGhostTimeCluster
    : public time_stepping::GhostTimeClusterWithCopy<Mpi::DataTransferMode::CopyInCopyOutHost> {
  public:
  using GhostTimeClusterWithCopy::GhostTimeClusterWithCopy;
  using GhostTimeClusterWithCopy::receiveGhostLayer;
  using GhostTimeClusterWithCopy::sendCopyLayer;
  using GhostTimeClusterWithCopy::testForCopyLayerSends;
  using GhostTimeClusterWithCopy::testForGhostLayerReceives;
};

// what rank `sender` sends in `region` at `step`; exact in double precision, and never zero (as
// fresh pinned memory is) or -1 (as the ghost regions are before an exchange)
inline double haloValue(int sender, std::size_t region, int step, std::size_t index) {
  return static_cast<double>(sender + 1) * 1.0e7 + static_cast<double>(region) * 1.0e6 +
         static_cast<double>(step) * 1.0e5 + static_cast<double>(index);
}

inline void exchangeRing(bool persistent) {
  auto& device = ::device::DeviceInstance::getInstance();
  const int rank = seissol::Mpi::mpi.rank();
  const int size = seissol::Mpi::mpi.size();

  // region 0 goes to the right neighbor and region 1 to the left one (for up to two ranks, that is
  // the same rank, and only the tag tells the regions apart); region 1 is large, so that MPI
  // implementations usually do not send it eagerly
  constexpr std::size_t RegionCount = 2;
  constexpr int StepCount = 3;
  const std::array<std::size_t, RegionCount> regionSizes{1000, 40000};
  const std::array<int, RegionCount> sendTo{(rank + 1) % size, (rank + size - 1) % size};
  const std::array<int, RegionCount> receiveFrom{sendTo[1], sendTo[0]};

  std::array<void*, RegionCount> copyDevice{};
  std::array<void*, RegionCount> ghostDevice{};
  solver::HaloCommunication halo(1, std::vector<solver::RemoteClusterPair>(1));
  for (std::size_t region = 0; region < RegionCount; ++region) {
    const auto bytes = regionSizes[region] * sizeof(double);
    copyDevice[region] = device.api->allocGlobMem(bytes);
    ghostDevice[region] = device.api->allocGlobMem(bytes);
    halo[0][0].copy.emplace_back(
        copyDevice[region], regionSizes[region], RealType::F64, sendTo[region], region);
    halo[0][0].ghost.emplace_back(
        ghostDevice[region], regionSizes[region], RealType::F64, receiveFrom[region], region);
  }

  {
    HostGhostTimeCluster cluster(1.0, 1, 0, 0, "test", "test", halo, persistent);

    // several exchanges, so that the staging buffers need to be refreshed each time
    for (int step = 0; step < StepCount; ++step) {
      cluster.setTime(static_cast<double>(step));

      for (std::size_t region = 0; region < RegionCount; ++region) {
        const auto bytes = regionSizes[region] * sizeof(double);
        std::vector<double> copyHost(regionSizes[region]);
        for (std::size_t i = 0; i < copyHost.size(); ++i) {
          copyHost[i] = haloValue(rank, region, step, i);
        }
        const std::vector<double> ghostHost(regionSizes[region], -1.0);
        device.api->copyTo(copyDevice[region], copyHost.data(), bytes);
        device.api->copyTo(ghostDevice[region], ghostHost.data(), bytes);
      }
      device.api->syncDevice();

      cluster.receiveGhostLayer();
      cluster.sendCopyLayer();
      bool received = false;
      bool sent = false;
      while (!received || !sent) {
        received = cluster.testForGhostLayerReceives();
        sent = cluster.testForCopyLayerSends();
      }

      for (std::size_t region = 0; region < RegionCount; ++region) {
        std::vector<double> ghostHost(regionSizes[region]);
        device.api->copyFrom(
            ghostHost.data(), ghostDevice[region], regionSizes[region] * sizeof(double));
        std::size_t mismatches = 0;
        for (std::size_t i = 0; i < ghostHost.size(); ++i) {
          if (ghostHost[i] != haloValue(receiveFrom[region], region, step, i)) {
            ++mismatches;
          }
        }
        CAPTURE(step);
        CAPTURE(region);
        CHECK(mismatches == 0);
      }
    }

    cluster.finalize();
  }

  for (std::size_t region = 0; region < RegionCount; ++region) {
    device.api->freeGlobMem(copyDevice[region]);
    device.api->freeGlobMem(ghostDevice[region]);
  }
}

} // namespace ghosttimeclusterwithcopytest

using namespace ghosttimeclusterwithcopytest;

TEST_CASE("Ghost clusters exchange device halos through host memory" * doctest::test_suite("mpi")) {
  auto& device = ::device::DeviceInstance::getInstance();
  device.api->setDevice(0);
  device.api->initialize();

  for (const bool persistent : {true, false}) {
    CAPTURE(persistent);
    exchangeRing(persistent);
  }

  // as in main(); otherwise, the device is torn down only by static destructors at exit
  device.api->finalize();
}

} // namespace seissol::unit_test

#endif // ACL_DEVICE

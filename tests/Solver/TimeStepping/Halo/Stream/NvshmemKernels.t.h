// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#if defined(ACL_DEVICE) && defined(USE_SHMEM) && defined(SHMEM_NVSHMEM)

#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/Stream/NvshmemKernels.h"

#include <Device/device.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <nvshmem_host.h>

namespace seissol::unit_test {

TEST_CASE("Replayed recordings of NVSHMEM exchanges count on" * doctest::test_suite("solver")) {
  // One processing element that exchanges with itself, as the SHMEM exchange does with its peers:
  // count the exchange, signal the count, put data with a signal, and wait for both. A recording of
  // this has to count on in each replay (the host API of NVSHMEM would keep the recorded counts).
  auto& device = ::device::DeviceInstance::instance();
  // NVSHMEM needs the device to be set up, as SeisSol does before
  device.api().setDevice(0);
  device.api().initialize();
  REQUIRE(solver::nvshmem::initialize(Mpi::mpi.comm()) == 0);
  const int pe = nvshmem_my_pe();

  constexpr std::size_t Bytes = 1024;
  auto* window = nvshmem_malloc(Bytes);
  auto* signals = static_cast<std::uint64_t*>(nvshmem_malloc(2 * sizeof(std::uint64_t)));
  auto* clearToSend = &signals[0];
  auto* arrived = &signals[1];
  auto* count = static_cast<std::uint64_t*>(device.api().allocGlobMem(sizeof(std::uint64_t)));
  auto* data = device.api().allocGlobMem(Bytes);
  REQUIRE(window != nullptr);
  REQUIRE(signals != nullptr);
  const std::array<std::uint64_t, 2> zeros{0, 0};
  device.api().copyTo(signals, zeros.data(), sizeof(zeros));
  device.api().copyTo(count, zeros.data(), sizeof(std::uint64_t));
  nvshmem_barrier_all();

  void* stream = device.api().createStream();
  const auto exchange = [&]() {
    solver::nvshmem::advance(count, stream);
    solver::nvshmem::signalCount(clearToSend, count, pe, stream);
    solver::nvshmem::waitCount(clearToSend, count, 1, stream);
    solver::nvshmem::putSignal(window, data, Bytes, arrived, pe, stream);
    solver::nvshmem::quiet(stream);
    solver::nvshmem::waitCount(arrived, count, 1, stream);
    REQUIRE(solver::nvshmem::lastError() == nullptr);
  };

  // once as usual, once recorded (which does not run it), replayed, and once more as usual
  constexpr std::uint64_t Replays = 5;
  exchange();
  auto graph = device.api().streamBeginCapture({stream});
  REQUIRE(graph.isInitialized());
  exchange();
  device.api().streamEndCapture(graph);
  for (std::uint64_t i = 0; i < Replays; ++i) {
    device.api().launchGraph(graph, stream);
  }
  exchange();
  device.api().syncStreamWithHost(stream);

  std::array<std::uint64_t, 2> signalValues{};
  std::uint64_t countValue = 0;
  device.api().copyFrom(signalValues.data(), signals, sizeof(signalValues));
  device.api().copyFrom(&countValue, count, sizeof(countValue));
  CHECK(countValue == Replays + 2);
  CHECK(signalValues[0] == Replays + 2);
  CHECK(signalValues[1] == Replays + 2);

  graph.reset();
  device.api().destroyGenericStream(stream);
  device.api().freeGlobMem(data);
  device.api().freeGlobMem(count);
  nvshmem_free(signals);
  nvshmem_free(window);
  solver::nvshmem::finalize();
}

} // namespace seissol::unit_test

#endif

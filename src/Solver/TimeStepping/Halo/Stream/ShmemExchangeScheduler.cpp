// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ShmemExchangeScheduler.h"

#include "Solver/TimeStepping/Halo/Stream/ExchangeScheduler.h"
#include "Solver/TimeStepping/Halo/Stream/StreamExchangeScheduler.h"

#include <cstddef>
#include <cstdint>
#include <utils/logger.h>
#include <vector>

#if defined(ACL_DEVICE) && defined(USE_SHMEM)

#include "Kernels/Precision.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/Halo/HaloCommunication.h"

#include <Device/device.h>
#include <algorithm>
#include <mpi.h>

#ifdef SHMEM_NVSHMEM
#include <nvshmem_host.h>
#endif
#ifdef SHMEM_ROCSHMEM
#include <rocshmem/rocshmem.hpp>
#endif
#ifdef SHMEM_ISHMEM
#include <ishmem.h>
#include <ishmemx.h>
#include <sycl/sycl.hpp>
#endif

namespace seissol::time_stepping {

namespace {

// the operations the protocol needs, for each library

#ifdef SHMEM_NVSHMEM
void shmemInit(MPI_Comm comm) {
  nvshmemx_init_attr_t attr{};
  nvshmemx_set_attr_mpi_comm_args(&comm, &attr);
  nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
}
void shmemFinalize() { nvshmem_finalize(); }
int shmemPe() { return nvshmem_my_pe(); }
void* shmemMalloc(std::size_t bytes) { return nvshmem_malloc(bytes); }
void shmemFree(void* pointer) { nvshmem_free(pointer); }
void shmemBarrierAll() { nvshmem_barrier_all(); }
void signalOnStream(std::uint64_t* signal, std::uint64_t value, int pe, void* stream) {
  nvshmemx_signal_op_on_stream(
      signal, value, NVSHMEM_SIGNAL_SET, pe, static_cast<cudaStream_t>(stream));
}
void waitOnStream(std::uint64_t* signal, std::uint64_t value, void* stream) {
  nvshmemx_signal_wait_until_on_stream(
      signal, NVSHMEM_CMP_GE, value, static_cast<cudaStream_t>(stream));
}
void putSignalOnStream(void* destination,
                       const void* source,
                       std::size_t bytes,
                       std::uint64_t* signal,
                       int pe,
                       void* stream) {
  nvshmemx_putmem_signal_nbi_on_stream(destination,
                                       source,
                                       bytes,
                                       signal,
                                       1,
                                       NVSHMEM_SIGNAL_ADD,
                                       pe,
                                       static_cast<cudaStream_t>(stream));
}
void quietOnStream(void* stream) { nvshmemx_quiet_on_stream(static_cast<cudaStream_t>(stream)); }
#endif

#ifdef SHMEM_ROCSHMEM
void shmemInit(MPI_Comm comm) {
  rocshmem::rocshmem_init_attr_t attr{};
  rocshmem::rocshmem_set_attr_mpi_comm_args(&comm, &attr);
  rocshmem::rocshmem_init_attr(rocshmem::ROCSHMEM_INIT_WITH_MPI_COMM, &attr);
}
void shmemFinalize() { rocshmem::rocshmem_finalize(); }
int shmemPe() { return rocshmem::rocshmem_my_pe(); }
void* shmemMalloc(std::size_t bytes) { return rocshmem::rocshmem_malloc(bytes); }
void shmemFree(void* pointer) { rocshmem::rocshmem_free(pointer); }
void shmemBarrierAll() { rocshmem::rocshmem_barrier_all(); }
void signalOnStream(std::uint64_t* signal, std::uint64_t value, int pe, void* stream) {
  rocshmem::rocshmem_signal_op_on_stream(
      signal, value, rocshmem::ROCSHMEM_SIGNAL_SET, pe, static_cast<hipStream_t>(stream));
}
void waitOnStream(std::uint64_t* signal, std::uint64_t value, void* stream) {
  rocshmem::rocshmem_signal_wait_until_on_stream(
      signal, rocshmem::ROCSHMEM_CMP_GE, value, static_cast<hipStream_t>(stream));
}
void putSignalOnStream(void* destination,
                       const void* source,
                       std::size_t bytes,
                       std::uint64_t* signal,
                       int pe,
                       void* stream) {
  rocshmem::rocshmem_putmem_signal_nbi_on_stream(destination,
                                                 source,
                                                 bytes,
                                                 signal,
                                                 1,
                                                 rocshmem::ROCSHMEM_SIGNAL_ADD,
                                                 pe,
                                                 static_cast<hipStream_t>(stream));
}
void quietOnStream(void* stream) {
  rocshmem::rocshmem_quiet_on_stream(static_cast<hipStream_t>(stream));
}
#endif

#ifdef SHMEM_ISHMEM
sycl::queue& queue(void* stream) { return *static_cast<sycl::queue*>(stream); }
void shmemInit(MPI_Comm comm) {
  ishmemx_attr_t attr{};
  attr.runtime = ISHMEMX_RUNTIME_MPI;
  attr.mpi_comm = &comm;
  ishmemx_init_attr(&attr);
}
void shmemFinalize() { ishmem_finalize(); }
int shmemPe() { return ishmem_my_pe(); }
void* shmemMalloc(std::size_t bytes) { return ishmem_malloc(bytes); }
void shmemFree(void* pointer) { ishmem_free(pointer); }
void shmemBarrierAll() { ishmem_barrier_all(); }
void signalOnStream(std::uint64_t* signal, std::uint64_t value, int pe, void* stream) {
  ishmemx_signal_op_on_queue(signal, value, ISHMEM_SIGNAL_SET, pe, queue(stream));
}
void waitOnStream(std::uint64_t* signal, std::uint64_t value, void* stream) {
  ishmemx_signal_wait_until_on_queue(signal, ISHMEM_CMP_GE, value, queue(stream));
}
void putSignalOnStream(void* destination,
                       const void* source,
                       std::size_t bytes,
                       std::uint64_t* signal,
                       int pe,
                       void* stream) {
  ishmemx_putmem_signal_nbi_on_queue(
      destination, source, bytes, signal, 1, ISHMEM_SIGNAL_ADD, pe, queue(stream));
}
void quietOnStream(void* stream) { ishmemx_quiet_on_queue(queue(stream)); }
#endif

device::DeviceInstance& deviceInstance() { return device::DeviceInstance::getInstance(); }

std::size_t bytesOf(const solver::RemoteCluster& region) {
  return region.size * sizeOfRealType(region.datatype);
}

// the peers of the regions, each once, in the order of the regions
std::vector<int> peersOf(const std::vector<solver::RemoteCluster>& regions) {
  std::vector<int> peers;
  for (const auto& region : regions) {
    if (std::find(peers.begin(), peers.end(), region.rank) == peers.end()) {
      peers.push_back(region.rank);
    }
  }
  return peers;
}

constexpr std::size_t Alignment = 256;

} // namespace

ShmemExchangeScheduler::ShmemExchangeScheduler(std::size_t clusterCount)
    : StreamExchangeScheduler(clusterCount, LaunchOrder::Global), rank_(Mpi::mpi.rank()),
      size_(Mpi::mpi.size()) {
  shmemInit(Mpi::mpi.comm());
  if (shmemPe() != rank_) {
    logError() << "The SHMEM processing elements are not numbered as the MPI ranks.";
  }
}

ShmemExchangeScheduler::~ShmemExchangeScheduler() {
  synchronize();
  shmemBarrierAll();
  if (window_ != nullptr) {
    shmemFree(arrived_);
    shmemFree(clearToSend_);
    shmemFree(window_);
  }
  shmemFinalize();
}

std::size_t ShmemExchangeScheduler::signalIndex(std::size_t from, std::size_t to, int peer) const {
  return (from * clusterCount() + to) * static_cast<std::size_t>(size_) +
         static_cast<std::size_t>(peer);
}

void ShmemExchangeScheduler::added(const ScheduledTransport& transport) {
  transports_.push_back(&transport);
  auto& offsets = ghostOffsets_[&transport];
  for (const auto& region : transport.regions().ghost) {
    offsets.push_back(windowSize_);
    windowSize_ += (bytesOf(region) + Alignment - 1) / Alignment * Alignment;
  }
}

void ShmemExchangeScheduler::prepare() {
  // symmetric memory has the same size on all processes
  std::uint64_t windowSize = std::max<std::uint64_t>(windowSize_, Alignment);
  MPI_Allreduce(MPI_IN_PLACE, &windowSize, 1, MPI_UINT64_T, MPI_MAX, Mpi::mpi.comm());
  window_ = static_cast<char*>(shmemMalloc(windowSize));

  const auto signalCount = clusterCount() * clusterCount() * static_cast<std::size_t>(size_);
  clearToSend_ = static_cast<std::uint64_t*>(shmemMalloc(signalCount * sizeof(std::uint64_t)));
  arrived_ = static_cast<std::uint64_t*>(shmemMalloc(signalCount * sizeof(std::uint64_t)));
  if (window_ == nullptr || clearToSend_ == nullptr || arrived_ == nullptr) {
    logError() << "Could not allocate the symmetric memory for the halo exchange.";
  }
  const std::vector<std::uint64_t> zeros(signalCount, 0);
  deviceInstance().api->copyTo(clearToSend_, zeros.data(), signalCount * sizeof(std::uint64_t));
  deviceInstance().api->copyTo(arrived_, zeros.data(), signalCount * sizeof(std::uint64_t));
  shmemBarrierAll();

  // the receivers tell the senders where to put their data, matched as the data would be
  std::vector<MPI_Request> requests;
  for (const auto* transport : transports_) {
    const auto& ghostRegions = transport->regions().ghost;
    auto& offsets = ghostOffsets_.at(transport);
    for (std::size_t i = 0; i < ghostRegions.size(); ++i) {
      requests.emplace_back();
      MPI_Isend(&offsets[i],
                1,
                MPI_UINT64_T,
                ghostRegions[i].rank,
                static_cast<int>(ghostRegions[i].tag),
                Mpi::mpi.comm(),
                &requests.back());
    }
    const auto& copyRegions = transport->regions().copy;
    auto& remote = remoteOffsets_[transport];
    remote.resize(copyRegions.size());
    for (std::size_t i = 0; i < copyRegions.size(); ++i) {
      requests.emplace_back();
      MPI_Irecv(&remote[i],
                1,
                MPI_UINT64_T,
                copyRegions[i].rank,
                static_cast<int>(copyRegions[i].tag),
                Mpi::mpi.comm(),
                &requests.back());
    }
  }
  MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
}

void ShmemExchangeScheduler::enqueueGroup(std::size_t slot,
                                          std::size_t from,
                                          std::size_t to,
                                          std::size_t exchange,
                                          const ScheduledTransport* sender,
                                          const ScheduledTransport* receiver) {
  void* current = stream(slot);
  const std::uint64_t count = exchange + 1;

  // the staging window is free again: the group before has copied its data out
  if (receiver != nullptr) {
    for (const auto peer : peersOf(receiver->regions().ghost)) {
      signalOnStream(&clearToSend_[signalIndex(from, to, rank_)], count, peer, current);
    }
  }

  if (sender != nullptr) {
    const auto& regions = sender->regions().copy;
    const auto& remote = remoteOffsets_.at(sender);
    for (const auto peer : peersOf(regions)) {
      waitOnStream(&clearToSend_[signalIndex(from, to, peer)], count, current);
      for (std::size_t i = 0; i < regions.size(); ++i) {
        if (regions[i].rank == peer) {
          putSignalOnStream(window_ + remote[i],
                            regions[i].data,
                            bytesOf(regions[i]),
                            &arrived_[signalIndex(from, to, rank_)],
                            peer,
                            current);
        }
      }
    }
    // the copy layer may overwrite its data once the puts have left
    quietOnStream(current);
  }

  if (receiver != nullptr) {
    const auto& regions = receiver->regions().ghost;
    const auto& offsets = ghostOffsets_.at(receiver);
    for (const auto peer : peersOf(regions)) {
      const auto regionsFromPeer = static_cast<std::uint64_t>(std::count_if(
          regions.begin(), regions.end(), [&](const auto& region) { return region.rank == peer; }));
      waitOnStream(&arrived_[signalIndex(from, to, peer)], count * regionsFromPeer, current);
    }
    for (std::size_t i = 0; i < regions.size(); ++i) {
      deviceInstance().api->copyBetweenAsync(
          regions[i].data, window_ + offsets[i], bytesOf(regions[i]), current);
    }
  }
}

} // namespace seissol::time_stepping

#else

namespace seissol::time_stepping {

ShmemExchangeScheduler::ShmemExchangeScheduler(std::size_t clusterCount)
    : StreamExchangeScheduler(clusterCount, LaunchOrder::Global) {
  logError() << "This build of SeisSol does not support exchanging the halo data with SHMEM.";
}

ShmemExchangeScheduler::~ShmemExchangeScheduler() = default;

void ShmemExchangeScheduler::prepare() {}

void ShmemExchangeScheduler::added(const ScheduledTransport& /*transport*/) {}

void ShmemExchangeScheduler::enqueueGroup(std::size_t /*slot*/,
                                          std::size_t /*from*/,
                                          std::size_t /*to*/,
                                          std::size_t /*exchange*/,
                                          const ScheduledTransport* /*sender*/,
                                          const ScheduledTransport* /*receiver*/) {}

std::size_t ShmemExchangeScheduler::signalIndex(std::size_t /*from*/,
                                                std::size_t /*to*/,
                                                int /*peer*/) const {
  return 0;
}

} // namespace seissol::time_stepping

#endif

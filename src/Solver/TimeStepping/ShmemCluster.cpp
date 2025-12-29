// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ShmemCluster.h"

#include "Memory/MemoryAllocator.h"
#include "Parallel/MPI.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <memory>
#include <mpi.h>

#ifdef USE_SHMEM

#ifdef USE_NVSHMEM
#include <nvshmem.h>
#include <nvshmemx.h>

using StreamT = cudaStream_t;
#endif

#ifdef USE_ROCSHMEM
#include <rocshmem.h>

using StreamT = hipStream_t;
#endif

#ifdef USE_ISHMEM
#include <ishmem.h>
#include <ishmemx.h>
#include <sycl/sycl.hpp>
using StreamT = sycl::queue*;
#endif

namespace seissol {
namespace {

void putOnStream(void* remote,
                 const void* local,
                 std::size_t size,
                 int rank,
                 RealType datatype,
                 uint64_t* signal,
                 void* stream) {
  StreamT nativeStream = static_cast<StreamT>(stream);
#ifdef USE_NVSHMEM
  if (datatype == RealType::F64) {
    nvhsmemx_double_put_signal_nbi_on_stream(
        remote, local, size, signal, 1, NVSHMEM_SIGNAL_ADD, rank, nativeStream);
  } else if (datatype == RealType::F32) {
    nvhsmemx_float_put_signal_nbi_on_stream(
        remote, local, size, signal, 1, NVSHMEM_SIGNAL_ADD, rank, nativeStream);
  }
#endif
}

void synchronize(void* stream) {
  StreamT nativeStream = static_cast<StreamT>(stream);
#ifdef USE_NVSHMEM
  nvshmemx_quiet_on_stream(nativeStream);
#endif
}

} // namespace
} // namespace seissol

#else

namespace seissol {
namespace {

void putOnStream(void* remote,
                 const void* local,
                 std::size_t size,
                 int rank,
                 RealType datatype,
                 uint64_t* signal,
                 void* stream) {}

void synchronize(void* stream) {}

} // namespace
} // namespace seissol

#endif

namespace seissol::time_stepping {

void ShmemCluster::launchSend() {
  for (std::size_t i = 0; i < meshStructure.copy.size(); ++i) {
    const auto& remote = meshStructure.copy[i];
    putOnStream(remotePointers[i],
                remote.data,
                remote.size,
                remote.rank,
                remote.datatype,
                signal.data(),
                stream->stream());
  }
}

bool ShmemCluster::mayCorrect() {
  synchronize(stream->stream());
  stream->wait();
  return AbstractTimeCluster::mayCorrect() && signal[0] == meshStructure.ghost.size();
}

void ShmemCluster::handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) {
  launchSend();
}

void ShmemCluster::handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) {}

void ShmemCluster::printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) {
  logError() << "Shmem timeout. More info TBD.";
}

void ShmemCluster::start() {}

void ShmemCluster::predict() {}
void ShmemCluster::correct() {}

ShmemCluster::ShmemCluster(double maxTimeStepSize,
                           int timeStepRate,
                           int globalTimeClusterId,
                           int otherGlobalTimeClusterId,
                           const seissol::solver::HaloCommunication& meshStructure,
                           bool /*persistent*/)
    : AbstractTimeCluster(
          maxTimeStepSize, timeStepRate, isDeviceOn() ? Executor::Device : Executor::Host),
      globalClusterId(globalTimeClusterId), otherGlobalClusterId(otherGlobalTimeClusterId),
      meshStructure(meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId)),
      signal(1, memory::Memkind::Shmem) {

  const auto& local = meshStructure.at(globalTimeClusterId).at(otherGlobalTimeClusterId);

  remotePointers.resize(local.ghost.size());
  std::vector<MPI_Request> requests(remotePointers.size());
  for (std::size_t i = 0; i < remotePointers.size(); ++i) {
    const auto& cluster = local.ghost[i];

    // send pointer as uintptr_t
    MPI_Isendrecv(static_cast<const void*>(&cluster.data),
                  1,
                  Mpi::castToMpiType<std::uintptr_t>(),
                  cluster.rank,
                  cluster.tag,
                  static_cast<void*>(&remotePointers[i]),
                  1,
                  Mpi::castToMpiType<std::uintptr_t>(),
                  cluster.rank,
                  cluster.tag,
                  Mpi::mpi.comm(),
                  &requests[i]);
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  signal[0] = 0;
}

ShmemCluster::~ShmemCluster() = default;

std::string ShmemCluster::description() const {
  return "comm-" + std::to_string(globalClusterId) + "-" + std::to_string(otherGlobalClusterId);
}

} // namespace seissol::time_stepping

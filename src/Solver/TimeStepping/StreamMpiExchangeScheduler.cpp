// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "StreamMpiExchangeScheduler.h"

#include "Solver/TimeStepping/ExchangeScheduler.h"
#include "Solver/TimeStepping/StreamExchangeScheduler.h"

#include <cstddef>
#include <memory>
#include <utils/logger.h>
#include <vector>

#if defined(ACL_DEVICE) && defined(USE_STREAM_MPI)

#include "Common/Typedefs.h"
#include "Kernels/Common.h"
#include "Parallel/MPI.h"
#include "Solver/TimeStepping/HaloCommunication.h"

#include <mpi.h>

#if !defined(STREAM_MPI_MPICH) && !defined(STREAM_MPI_CRAY)
#error "Stream-aware MPI needs STREAM_MPI_MPICH or STREAM_MPI_CRAY."
#endif

namespace seissol::time_stepping {

namespace {

void check(int result, const char* call) {
  if (result != MPI_SUCCESS) {
    logError() << "The stream-aware MPI call" << call << "failed with the error code" << result;
  }
}

#ifdef STREAM_MPI_MPICH
// the name MPICH knows the stream type of the device backend by
const char* streamType() {
  if constexpr (Backend == DeviceBackend::Cuda) {
    return "cudaStream_t";
  } else if constexpr (Backend == DeviceBackend::Hip) {
    return "hipStream_t";
  } else {
    logError() << "MPICH enqueues operations only on CUDA and HIP streams.";
    return "";
  }
}
#endif

} // namespace

struct StreamMpiExchangeScheduler::Slot {
#ifdef STREAM_MPI_MPICH
  MPIX_Stream stream{MPIX_STREAM_NULL};
  MPI_Comm comm{MPI_COMM_NULL};
#endif
#ifdef STREAM_MPI_CRAY
  MPIX_Queue queue{};
#endif
};

StreamMpiExchangeScheduler::StreamMpiExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : StreamExchangeScheduler(clusterCount, order) {
  slots_.resize(order == LaunchOrder::Global ? 1 : clusterCount * clusterCount);
  for (const auto current : usedSlots()) {
    auto slotData = std::make_unique<Slot>();
    void* deviceStream = stream(current);
#ifdef STREAM_MPI_MPICH
    MPI_Info info{};
    check(MPI_Info_create(&info), "MPI_Info_create");
    check(MPI_Info_set(info, "type", streamType()), "MPI_Info_set");
    check(MPIX_Info_set_hex(info, "value", &deviceStream, sizeof(deviceStream)),
          "MPIX_Info_set_hex");
    check(MPIX_Stream_create(info, &slotData->stream), "MPIX_Stream_create");
    check(MPI_Info_free(&info), "MPI_Info_free");
    check(MPIX_Stream_comm_create(Mpi::mpi.comm(), slotData->stream, &slotData->comm),
          "MPIX_Stream_comm_create");
#endif
#ifdef STREAM_MPI_CRAY
    check(MPIX_Create_queue(Mpi::mpi.comm(), deviceStream, &slotData->queue), "MPIX_Create_queue");
#endif
    slots_[current] = std::move(slotData);
  }
}

StreamMpiExchangeScheduler::~StreamMpiExchangeScheduler() {
  synchronize();
  for (auto& slotData : slots_) {
    if (slotData != nullptr) {
#ifdef STREAM_MPI_MPICH
      check(MPI_Comm_free(&slotData->comm), "MPI_Comm_free");
      check(MPIX_Stream_free(&slotData->stream), "MPIX_Stream_free");
#endif
#ifdef STREAM_MPI_CRAY
      check(MPIX_Free_queue(slotData->queue), "MPIX_Free_queue");
#endif
    }
  }
}

void StreamMpiExchangeScheduler::enqueueGroup(std::size_t slot,
                                              std::size_t /*from*/,
                                              std::size_t /*to*/,
                                              std::size_t /*exchange*/,
                                              const ScheduledTransport* sender,
                                              const ScheduledTransport* receiver) {
  auto& slotData = *slots_[slot];
  std::vector<MPI_Request> requests;
  if (sender != nullptr) {
    for (const auto& region : sender->regions().copy) {
      MPI_Request request{MPI_REQUEST_NULL};
#ifdef STREAM_MPI_MPICH
      check(MPIX_Isend_enqueue(region.data,
                               static_cast<int>(region.size),
                               Mpi::precisionToMpiType(region.datatype),
                               region.rank,
                               static_cast<int>(region.tag),
                               slotData.comm,
                               &request),
            "MPIX_Isend_enqueue");
#endif
#ifdef STREAM_MPI_CRAY
      check(MPIX_Enqueue_send(region.data,
                              static_cast<int>(region.size),
                              Mpi::precisionToMpiType(region.datatype),
                              region.rank,
                              static_cast<int>(region.tag),
                              slotData.queue,
                              &request),
            "MPIX_Enqueue_send");
#endif
      requests.push_back(request);
    }
  }
  if (receiver != nullptr) {
    for (const auto& region : receiver->regions().ghost) {
      MPI_Request request{MPI_REQUEST_NULL};
#ifdef STREAM_MPI_MPICH
      check(MPIX_Irecv_enqueue(region.data,
                               static_cast<int>(region.size),
                               Mpi::precisionToMpiType(region.datatype),
                               region.rank,
                               static_cast<int>(region.tag),
                               slotData.comm,
                               &request),
            "MPIX_Irecv_enqueue");
#endif
#ifdef STREAM_MPI_CRAY
      check(MPIX_Enqueue_recv(region.data,
                              static_cast<int>(region.size),
                              Mpi::precisionToMpiType(region.datatype),
                              region.rank,
                              static_cast<int>(region.tag),
                              slotData.queue,
                              &request),
            "MPIX_Enqueue_recv");
#endif
      requests.push_back(request);
    }
  }
  if (requests.empty()) {
    return;
  }
  // the stream continues once all operations of the group have completed
#ifdef STREAM_MPI_MPICH
  check(
      MPIX_Waitall_enqueue(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
      "MPIX_Waitall_enqueue");
#endif
#ifdef STREAM_MPI_CRAY
  check(MPIX_Enqueue_start(slotData.queue), "MPIX_Enqueue_start");
  check(MPIX_Enqueue_wait(slotData.queue), "MPIX_Enqueue_wait");
#endif
}

} // namespace seissol::time_stepping

#else

namespace seissol::time_stepping {

struct StreamMpiExchangeScheduler::Slot {};

StreamMpiExchangeScheduler::StreamMpiExchangeScheduler(std::size_t clusterCount, LaunchOrder order)
    : StreamExchangeScheduler(clusterCount, order) {
  logError() << "This build of SeisSol does not support exchanging the halo data with"
             << "stream-aware MPI.";
}

StreamMpiExchangeScheduler::~StreamMpiExchangeScheduler() = default;

void StreamMpiExchangeScheduler::enqueueGroup(std::size_t /*slot*/,
                                              std::size_t /*from*/,
                                              std::size_t /*to*/,
                                              std::size_t /*exchange*/,
                                              const ScheduledTransport* /*sender*/,
                                              const ScheduledTransport* /*receiver*/) {}

} // namespace seissol::time_stepping

#endif

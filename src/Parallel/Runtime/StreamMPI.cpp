// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "StreamMPI.h"

#include <mpi.h>

/* ---------------------------------------------------------------------
 * Feature detection
 * -------------------------------------------------------------------*/
#if defined(CRAY_MPICH_VERSION)
/* HPE Cray MPICH: MPIX_Enqueue_* API */
#define STREAM_MPI_BACKEND_CRAY 1
#elif defined(MPICH) && defined(MPICH_NUMVERSION) && (MPICH_NUMVERSION >= 40100000)
/* Argonne MPICH 4.1+: MPIX_*_enqueue API */
#define STREAM_MPI_BACKEND_MPICH 1
#else
/* Fallback: classic CUDA-aware MPI with stream sync */
#define STREAM_MPI_BACKEND_FALLBACK 1
#endif

namespace seissol {

struct MPIStreamData {
#if defined(STREAM_MPI_BACKEND_CRAY)
  MPIX_Queue queue;
#elif defined(STREAM_MPI_BACKEND_MPICH)
  MPIX_Stream stream;
  MPI_Comm scomm;
#endif
};

MPIStream::MPIStream(void* stream) : stream_(stream) { init(); }

void MPIStream::init() {
#if defined(STREAM_MPI_BACKEND_CRAY)
  MPIX_Create_queue(comm_, stream_, &data_->queue);
#elif defined(STREAM_MPI_BACKEND_MPICH)
  MPI_Info info;
  MPI_Info_create(&info);
  /* the stream as hex info */
  MPIX_Info_set_hex(info, "type", "cudaStream_t", sizeof("cudaStream_t") - 1);
  MPIX_Info_set_hex(info, "value", (void*)&stream_, sizeof(s));
  MPIX_Stream_create(info, &data_->stream);
  MPI_Info_free(&info);
  MPIX_Stream_comm_create(comm_, data_->stream, &data_->scomm);
#endif
}

void MPIStream::send(
    void* data, std::size_t count, MPI_Datatype datatype, int other, int tag, MPI_Request request) {
#if defined(STREAM_MPI_BACKEND_CRAY)
  MPIX_Enqueue_send(data, count, datatype, other, tag, data_->queue, request);
#elif defined(STREAM_MPI_BACKEND_MPICH)
  MPIX_Isend_enqueue(data, count, datatype, other, tag, data_->scomm, request);
#endif
  requests_->push_back(request);
}

void MPIStream::recv(
    void* data, std::size_t count, MPI_Datatype datatype, int other, int tag, MPI_Request request) {
#if defined(STREAM_MPI_BACKEND_CRAY)
  MPIX_Enqueue_recv(data, count, datatype, other, tag, data_->queue, request);
#elif defined(STREAM_MPI_BACKEND_MPICH)
  MPIX_Irecv_enqueue(data, count, datatype, other, tag, data_->scomm, request);
#endif
  requests_->push_back(request);
}

void MPIStream::commit() {
#if defined(STREAM_MPI_BACKEND_CRAY)
  MPIX_Enqueue_start(data_->queue);
#endif
}

void MPIStream::complete() {
#if defined(STREAM_MPI_BACKEND_CRAY)
  MPIX_Enqueue_wait(data_->queue);
#elif defined(STREAM_MPI_BACKEND_MPICH)
  MPIX_Waitall_enqueue(requests_->size(), requests_->data(), data_->stream);
#endif
  requests_->clear();
}

} // namespace seissol

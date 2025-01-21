// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Stream.h"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.h"

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
#ifdef SEISSOL_KERNELS_SYCL
  device().api->recordEventOnStream(forkEventSycl, streamPtr);
  device().api->syncStreamWithEvent(queuePtr, forkEventSycl);
#else
  sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
  auto* localForkEventSycl{forkEventSycl};
  device().api->recordEventOnStream(localForkEventSycl, streamPtr);
  syclNativeOperation(*queue, true, [=](void* stream) {
    device().api->syncStreamWithEvent(stream, localForkEventSycl);
  });
#endif
}
void StreamRuntime::syncFromSycl(void* queuePtr) {
#ifdef SEISSOL_KERNELS_SYCL
  device().api->recordEventOnStream(joinEventSycl, queuePtr);
  device().api->syncStreamWithEvent(streamPtr, joinEventSycl);
#else
  sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
  auto syclEvent = syclNativeOperation(*queue, true, [=](void* stream) {
    device().api->recordEventOnStream(joinEventSycl, stream);
  });

  // needs a submission barrier here
  // a bit hacky right now; but it works
#if defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(SYCL_EXT_ACPP_ENQUEUE_CUSTOM_OPERATION)
  if (queue->get_context().AdaptiveCpp_runtime() != nullptr) {
    queue->get_context().AdaptiveCpp_runtime()->dag().flush_sync();
  }
#elif defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION)
  if (queue->get_context().hipSYCL_runtime() != nullptr) {
    queue->get_context().hipSYCL_runtime()->dag().flush_sync();
  }
#else
  // note that polling on the info value for "is not pending" may not advance the DAG
  // (at least in the case of AdaptiveCpp)
  syclEvent.wait();
#endif
  device().api->syncStreamWithEvent(streamPtr, joinEventSycl);
#endif

  /*
  // the following will not work, because SYCL may decide to postpone execution a bit
  // effectively, there may be no solution but removing SYCL entirely, or making SYCL the
  over-arching runtime
  // or, potentially OpenMP could play as glue, once it gains its interop functionality
  auto* localJoinEventSycl{joinEventSycl};
  syclNativeOperation(*queue, true, [=](void* stream) {
    device().api->recordEventOnStream(localJoinEventSycl, stream);
  });
  device().api->syncStreamWithEvent(streamPtr, localJoinEventSycl);
  */
}

} // namespace seissol::parallel::runtime

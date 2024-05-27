#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"

#include <semaphore.h>

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
  auto* localJoinEventSycl{joinEventSycl};
  auto event = syclNativeOperation(*queue, true, [=](void* stream) {
    device().api->recordEventOnStream(localJoinEventSycl, stream);
  });
  // needs a submission barrier here
#if defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION)
  // a bit hacky right now; but it works
  if (queue->get_context().hipSYCL_runtime() != nullptr) {
    queue->get_context().hipSYCL_runtime()->dag().flush_sync();
  }
#else
  // note that polling on the info value for "is not pending" may not advance the DAG
  // (at least in the case of AdaptiveCpp)
  event.wait();
#endif
  device().api->syncStreamWithEvent(streamPtr, localJoinEventSycl);
#endif
}

} // namespace seissol::parallel::runtime

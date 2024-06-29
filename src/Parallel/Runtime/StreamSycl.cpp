#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"

#include <semaphore.h>

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
  void* event = getEvent();
#ifdef SEISSOL_KERNELS_SYCL
  device().api->recordEventOnStream(event, streamPtr);
  device().api->syncStreamWithEvent(queuePtr, event);
#else
  sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
  device().api->recordEventOnStream(event, streamPtr);
  syclNativeOperation(
      *queue, true, [=](void* stream) { device().api->syncStreamWithEvent(stream, event); });
#endif
}
void StreamRuntime::syncFromSycl(void* queuePtr) {
  void* event = getEvent();
#ifdef SEISSOL_KERNELS_SYCL
  device().api->recordEventOnStream(event, queuePtr);
  device().api->syncStreamWithEvent(streamPtr, event);
#else
  sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
  auto syclEvent = syclNativeOperation(
      *queue, true, [=](void* stream) { device().api->recordEventOnStream(event, stream); });

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
  device().api->syncStreamWithEvent(streamPtr, event);
#endif
}

} // namespace seissol::parallel::runtime

#include "Stream.h"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.h"

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
  void* event = getEvent();
#ifdef SEISSOL_KERNELS_SYCL
  device().api->recordEventOnStream(event, streamPtr);
  device().api->syncStreamWithEvent(queuePtr, event);
#else
  sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
  auto* localForkEventSycl{event};
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
  queue->wait();
#endif
}

} // namespace seissol::parallel::runtime

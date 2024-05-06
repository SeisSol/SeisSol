#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
#ifdef SEISSOL_KERNELS_SYCL
  device().api->recordEventOnStream(localForkEventSycl, streamPtr);
  device().api->syncStreamWithEvent(queuePtr, localForkEventSycl);
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
  sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
  queue->wait();

  /*
  // the following will not work, because SYCL may decide to postpone execution a bit
  auto* localJoinEventSycl{joinEventSycl};
  syclNativeOperation(*queue, true, [=](void* stream) {
    device().api->recordEventOnStream(localJoinEventSycl, stream);
  });
  device().api->syncStreamWithEvent(streamPtr, localJoinEventSycl);
  */
}

} // namespace seissol::parallel::runtime

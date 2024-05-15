#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"

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
  queue->wait();
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

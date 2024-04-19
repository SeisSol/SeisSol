#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
    sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
    auto* localForkEventSycl{forkEventSycl};
    device().api->recordEventOnStream(localForkEventSycl, streamPtr);
    syclNativeOperation(*queue, true, [=](void* stream) {
        device().api->syncStreamWithEvent(stream, localForkEventSycl);
    });
}
void StreamRuntime::syncFromSycl(void* queuePtr) {
    sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
    auto* localJoinEventSycl{joinEventSycl};
    syclNativeOperation(*queue, true, [=](void* stream) {
        device().api->recordEventOnStream(localJoinEventSycl, stream);
    });
    device().api->syncStreamWithEvent(streamPtr, localJoinEventSycl);
}

} // namespace seissol::parallel::runtime

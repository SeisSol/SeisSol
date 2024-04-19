#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
    sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
    device().api->recordEventOnStream(forkEventSycl, streamPtr);
    syclNativeOperation(*queue, true, [=](void* stream) {
        device().api->syncStreamWithEvent(stream, forkEventSycl);
    });
}
void StreamRuntime::syncFromSycl(void* queuePtr) {
    sycl::queue* queue = static_cast<sycl::queue*>(queuePtr);
    syclNativeOperation(*queue, true, [=](void* stream) {
        device().api->recordEventOnStream(joinEventSycl, stream);
    });
    device().api->syncStreamWithEvent(streamPtr, joinEventSycl);
}

} // namespace seissol::parallel::runtime

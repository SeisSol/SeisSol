#include "Stream.h"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.h"

namespace seissol::parallel::runtime {

void StreamRuntime::syncToSycl(void* queuePtr) {
}
void StreamRuntime::syncFromSycl(void* queuePtr) {
}

} // namespace seissol::parallel::runtime

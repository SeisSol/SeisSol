#include "Stream.hpp"

#include "Parallel/AcceleratorDevice.h"
#include "Parallel/SyclInterop.hpp"
#include <omp.h>

namespace seissol::parallel::runtime {

// TODO: maybe use OMP events? (i.e. ints that have only in/out dependencies)

/*
omp_interop_t interop;
#pragma omp interop depend(inout:queue[0]) nowait init(targetsync: interop)
void* stream = omp_get_interop_ptr(interop, omp_ipr_targetsync, nullptr);
device().api->syncStreamWithEvent(stream, forkEventSycl);
#pragma omp interop depend(inout:queue[0]) destroy(interop)
*/

#ifdef ACL_DEVICE
void StreamRuntime::syncToOMP(omp_depend_t& dep) {
  omp_event_handle_t handle;
#pragma omp task detach(handle) depend(depobj : dep)
  {}
  device().api->streamHostFunction(stream(), [=]() { omp_fulfill_event(handle); });
}
void StreamRuntime::syncFromOMP(omp_depend_t& dep) {
  auto* event = joinEventSycl;
#pragma omp task depend(depobj : dep)
  { device().api->recordEventOnHost(event); }
  device().api->syncStreamWithEvent(streamPtr, event);
}
#endif

} // namespace seissol::parallel::runtime

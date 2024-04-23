#include "Stream.hpp"

#include <omp.h>

#ifdef ACL_DEVICE
#include "device.h"
#endif

namespace seissol::parallel::runtime {

/*
// unused code:
omp_interop_t interop;
#pragma omp interop depend(inout:queue[0]) nowait init(targetsync: interop)
void* stream = omp_get_interop_ptr(interop, omp_ipr_targetsync, nullptr);
device().api->syncStreamWithEvent(stream, forkEventSycl);
#pragma omp interop depend(inout:queue[0]) destroy(interop)
*/

// currently completely disabled... Needs a modern compiler (OMP 5.0-capable)
#ifdef ACL_DEVICE
#if 0
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
#endif

} // namespace seissol::parallel::runtime

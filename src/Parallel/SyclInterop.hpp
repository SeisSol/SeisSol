#pragma once

#include "utils/logger.h"
#include <Parallel/AcceleratorDevice.h>
#include <utility>

#if defined(__HIPSYCL_ENABLE_CUDA_TARGET__) || defined(__ACPP_ENABLE_CUDA_TARGET__)
#define SEISSOL_SYCL_BACKEND_CUDA
#endif
#if defined(__HIPSYCL_ENABLE_HIP_TARGET__) || defined(__ACPP_ENABLE_HIP_TARGET__)
#define SEISSOL_SYCL_BACKEND_HIP
#endif
#if defined(__HIPSYCL_ENABLE_SPIRV_TARGET__) || defined(__ACPP_ENABLE_SPIRV_TARGET__)
#define SEISSOL_SYCL_BACKEND_SPIRV
#endif

namespace seissol::parallel {
// this file here does more or less try to interop to CUDA or HIP.
// and yes, it's more or less a preprocessor chaos
#if defined(SYCL_EXT_ONEAPI_ENQUEUE_BARRIER) || defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) ||            \
    defined(ACPP_EXT_QUEUE_WAIT_LIST)
#define SEISSOL_SYCL_HAS_BARRIER
#endif
#if defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION) || \
    defined(ONEAPI)
#define SEISSOL_SYCL_HAS_HOST_TASK
#endif

template <typename F>
sycl::event syclNativeOperation(sycl::queue& queue, bool blocking, F&& function) {
#ifndef SEISSOL_SYCL_HAS_HOST_TASK
  logError() << "Requested SYCL native interop operation, but that is not supported.";
#endif
#ifndef SEISSOL_SYCL_HAS_BARRIER
  if (blocking) {
    logError() << "Requested blocking SYCL operation, but no blocking supported.";
  }
#endif
#ifdef SYCL_EXT_ONEAPI_ENQUEUE_BARRIER
  if (blocking) {
    queue.ext_oneapi_submit_barrier();
  }
#endif
#if defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) || defined(ACPP_EXT_QUEUE_WAIT_LIST)
  auto waitList = queue.get_wait_list();
#endif
  // NOTE: if some thread adds something to the queue here, we may have a problem
  return queue.submit([&](sycl::handler& h) {
#if defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) || defined(ACPP_EXT_QUEUE_WAIT_LIST)
    // use "barrier", if needed
    if (blocking) {
      h.depends_on(waitList);
    }
#endif
#if defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION)
    h.hipSYCL_enqueue_custom_operation([=](sycl::interop_handle& handle) {
#ifdef SEISSOL_SYCL_BACKEND_CUDA
      if (queue.get_device().get_backend() == sycl::backend::cuda) {
        auto stream = handle.get_native_queue<sycl::backend::cuda>();
        std::invoke(function, stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
      if (queue.get_device().get_backend() == sycl::backend::hip) {
        auto stream = handle.get_native_queue<sycl::backend::hip>();
        std::invoke(function, stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_ZE
      if (queue.get_device().get_backend() == sycl::backend::ze) {
        auto stream = handle.get_native_queue<sycl::backend::ze>();
        std::invoke(function, stream);
        return;
      }
#endif
      { logError() << "Unknown backend" << static_cast<int>(queue.get_device().get_backend()); }
    });
#endif
#ifdef ONEAPI
    h.host_task([=](sycl::interop_handle& handle) {
#ifdef SEISSOL_SYCL_BACKEND_CUDA
      if (queue.get_device().get_backend() == sycl::backend::cuda) {
        auto stream = handle.get_native<sycl::backend::cuda, sycl::queue>();
        std::invoke(function, stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
      if (queue.get_device().get_backend() == sycl::backend::hip) {
        auto stream = handle.get_native<sycl::backend::hip, sycl::queue>();
        std::invoke(function, stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_ZE
      if (queue.get_device().get_backend() == sycl::backend::ze) {
        auto stream = handle.get_native<sycl::backend::ze, sycl::queue>();
        std::invoke(function, stream);
        return;
      }
#endif
    });
#endif
  });
}

sycl::event syclQueueSynchronize(sycl::queue& queueFrom, sycl::queue& queueTo) {
#ifdef SYCL_EXT_ONEAPI_ENQUEUE_BARRIER
  auto eventFrom = queueFrom.ext_oneapi_submit_barrier();
  auto eventTo = queueTo.ext_oneapi_submit_barrier({eventFrom});
  return eventTo;
#elif defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) || defined(ACPP_EXT_QUEUE_WAIT_LIST)
  auto waitList1 = queueFrom.get_wait_list();
  auto waitList2 = queueTo.get_wait_list();
  auto queueEvent = queueTo.submit([&](sycl::handler& h) {
    h.depends_on(waitList1);
    h.depends_on(waitList2);
#if defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION)
    h.hipSYCL_enqueue_custom_operation([=](auto&) {});
#else
    h.single_task([=](auto&) {});
#endif
  });
  return queueEvent;
#else
  logError() << "Queue synchronization is not implemented for this SYCL implementation.";
#endif
}

} // namespace seissol::parallel

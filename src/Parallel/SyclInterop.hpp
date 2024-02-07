#pragma once

#include "utils/logger.h"
#include <utility>
#include <Parallel/AcceleratorDevice.h>

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
void syclNativeOperation(sycl::queue& queue, bool blocking, F&& function) {
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
  queue.submit([&](sycl::handler& h) {
#if defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) || defined(ACPP_EXT_QUEUE_WAIT_LIST)
    // use "barrier", if needed
    if (blocking) {
      h.depends_on(queue.get_wait_list());
    }
#endif
#if defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION)
    h.hipSYCL_enqueue_custom_operation([&](sycl::interop_handle& handle) {
#ifdef SEISSOL_SYCL_BACKEND_CUDA
      if (queue.get_device().get_backend() == sycl::backend::cuda) {
        auto stream = handle.get_native_queue<sycl::backend::cuda>();
        std::invoke(std::forward<F>(function), stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
      if (queue.get_device().get_backend() == sycl::backend::hip) {
        auto stream = handle.get_native_queue<sycl::backend::hip>();
        std::invoke(std::forward<F>(function), stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_ZE
      if (queue.get_device().get_backend() == sycl::backend::ze) {
        auto stream = handle.get_native_queue<sycl::backend::ze>();
        std::invoke(std::forward<F>(function), stream);
        return;
      }
#endif
      { logError() << "Unknown backend" << static_cast<int>(queue.get_device().get_backend()); }
    });
#endif
#ifdef ONEAPI
    h.host_task([=](sycl::interop_handle& handle) {
      auto stream = handle.get_native_queue();
      std::invoke(std::forward<F>(function), stream);
    });
#endif
  });
}
} // namespace seissol::parallel

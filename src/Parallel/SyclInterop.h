// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PARALLEL_SYCLINTEROP_H_
#define SEISSOL_SRC_PARALLEL_SYCLINTEROP_H_

#include "utils/logger.h"
#include <Parallel/AcceleratorDevice.h>
#include <utility>

#if defined(__HIPSYCL_ENABLE_CUDA_TARGET__) || defined(__ACPP_ENABLE_CUDA_TARGET__) ||             \
    defined(SYCL_BACKEND_CUDA) || defined(SYCL_EXT_ONEAPI_BACKEND_CUDA)
#define SEISSOL_SYCL_BACKEND_CUDA
#endif
#if defined(__HIPSYCL_ENABLE_HIP_TARGET__) || defined(__ACPP_ENABLE_HIP_TARGET__) ||               \
    defined(SYCL_BACKEND_HIP) || defined(SYCL_EXT_ONEAPI_BACKEND_HIP)
#define SEISSOL_SYCL_BACKEND_HIP
#endif
#if defined(__HIPSYCL_ENABLE_SPIRV_TARGET__) || defined(__ACPP_ENABLE_SPIRV_TARGET__)
#define SEISSOL_SYCL_BACKEND_SPIRV
#endif

namespace seissol::parallel {
// this file here does more or less try to interop to CUDA or HIP.
// and yes, it's more or less a preprocessor chaos
#if defined(SYCL_EXT_ONEAPI_ENQUEUE_BARRIER) || defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) ||            \
    defined(ACPP_EXT_QUEUE_WAIT_LIST) || defined(SYCL_EXT_ACPP_QUEUE_WAIT_LIST)
#define SEISSOL_SYCL_HAS_BARRIER
#endif

#if defined(__HIPSYCL_ENABLE_CUDA_TARGET__) || defined(__ACPP_ENABLE_CUDA_TARGET__) ||             \
    defined(SYCL_BACKEND_CUDA)
constexpr auto SyclBackendCUDA = sycl::backend::cuda;
#elif defined(SYCL_EXT_ONEAPI_BACKEND_CUDA)
constexpr auto SyclBackendCUDA = sycl::backend::ext_oneapi_cuda;
#elif defined(SEISSOL_SYCL_BACKEND_CUDA)
// assert that it may never happen
#error "SYCL CUDA backend found; but no type for it."
#endif

#if defined(__HIPSYCL_ENABLE_HIP_TARGET__) || defined(__ACPP_ENABLE_HIP_TARGET__) ||               \
    defined(SYCL_BACKEND_HIP)
constexpr auto SyclBackendHIP = sycl::backend::hip;
#elif defined(SYCL_EXT_ONEAPI_BACKEND_HIP)
constexpr auto SyclBackendHIP = sycl::backend::ext_oneapi_hip;
#elif defined(SEISSOL_SYCL_BACKEND_HIP)
// assert that it may never happen
#error "SYCL HIP backend found; but no type for it."
#endif

template <typename F>
sycl::event syclNativeOperation(sycl::queue& queue, bool blocking, F&& function) {
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
#if defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) || defined(ACPP_EXT_QUEUE_WAIT_LIST) ||                   \
    defined(SYCL_EXT_ACPP_QUEUE_WAIT_LIST)
  auto waitList = queue.get_wait_list();
#endif
  // NOTE: if some thread adds something to the queue here, we may have a problem
  return queue.submit([&](sycl::handler& h) {
#if defined(HIPSYCL_EXT_QUEUE_WAIT_LIST) || defined(ACPP_EXT_QUEUE_WAIT_LIST) ||                   \
    defined(SYCL_EXT_ACPP_QUEUE_WAIT_LIST)
    // use "barrier", if needed
    if (blocking) {
      h.depends_on(waitList);
    }
#endif
#if defined(HIPSYCL_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION) || \
    defined(SYCL_EXT_ACPP_ENQUEUE_CUSTOM_OPERATION)
#if defined(ACPP_EXT_ENQUEUE_CUSTOM_OPERATION) || defined(SYCL_EXT_ACPP_ENQUEUE_CUSTOM_OPERATION)
    h.AdaptiveCpp_enqueue_custom_operation([=](sycl::interop_handle& handle) {
#else
    h.hipSYCL_enqueue_custom_operation([=](sycl::interop_handle& handle) {
#endif
#ifdef SEISSOL_SYCL_BACKEND_CUDA
      if (queue.get_device().get_backend() == SyclBackendCUDA) {
        auto stream = handle.get_native_queue<SyclBackendCUDA>();
        std::invoke(function, stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
      if (queue.get_device().get_backend() == SyclBackendHIP) {
        auto stream = handle.get_native_queue<SyclBackendHIP>();
        std::invoke(function, stream);
        return;
      }
#endif
      {
        logError() << "Unknown backend" << static_cast<int>(queue.get_device().get_backend());
      }
    });
#else
    // we cannot take the fast path; so just submit a host task instead
    h.host_task([=](sycl::interop_handle& handle) {
#ifdef SEISSOL_SYCL_BACKEND_CUDA
      if (queue.get_device().get_backend() == SyclBackendCUDA) {
        auto stream = handle.get_native<SyclBackendCUDA, sycl::queue>();
        std::invoke(function, stream);
        return;
      }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
      if (queue.get_device().get_backend() == SyclBackendHIP) {
        auto stream = handle.get_native<SyclBackendHIP, sycl::queue>();
        std::invoke(function, stream);
        return;
      }
#endif
      {
        logError() << "Unknown backend" << static_cast<int>(queue.get_device().get_backend());
      }
    });
#endif
  });
}
} // namespace seissol::parallel

#endif // SEISSOL_SRC_PARALLEL_SYCLINTEROP_H_

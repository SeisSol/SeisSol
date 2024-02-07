#pragma once

#include <Parallel/AcceleratorDevice.h>
#include <Parallel/SyclInterop.hpp>
#include <optional>

#ifndef __DPCPP_COMPILER
namespace sycl = cl::sycl;
#endif

namespace seissol::parallel {
class SyclNoGraph {
  public:
  using StreamT = sycl::queue;

  bool canCapture() const { return false; }

  void beginCapture(StreamT& queue) {}
  void endCapture(StreamT& queue) {}
  void runCapture(StreamT& queue) {}
};

template <typename NativeStreamT, typename NativeGraphT, typename NativeInstanceT>
class SyclNativeGraph {
  public:
  using StreamT = sycl::queue;

  bool canCapture() const;

  void beginCapture(StreamT& queue);
  void endCapture(StreamT& queue);
  void runCapture(StreamT& queue);

  private:
  NativeGraphT graph;
  std::optional<NativeInstanceT> instance;
};
#ifdef SYCL_EXT_ONEAPI_GRAPH
class SyclOneapiGraph {
  public:
  using StreamT = sycl::queue;

  bool canCapture() const;

  void beginCapture(StreamT& queue);
  void endCapture(StreamT& queue);
  void runCapture(StreamT& queue);

  private:
  sycl::ext::oneapi::experimental::command_graph graph;
  std::optional<sycl::ext::oneapi::experimental::command_graph<
      sycl::ext::oneapi::experimental::graph_state::executable>>
      instance;
};
using SyclGraph = SyclOneapiGraph;
#elif defined(SEISSOL_SYCL_BACKEND_CUDA)
using SyclGraph = SyclNativeGraph<cudaStream_t, cudaGraph_t, cudaGraphExec_t>;
#elif defined(SEISSOL_SYCL_BACKEND_HIP)
using SyclGraph = SyclNativeGraph<hipStream_t, hipGraph_t, hipGraphInstance_t>;
#else
using SyclGraph = SyclNoGraph;
#endif
} // namespace seissol::parallel

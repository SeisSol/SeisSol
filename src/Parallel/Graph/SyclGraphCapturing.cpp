#include "SyclGraphCapturing.hpp"

#include "Parallel/SyclInterop.hpp"

namespace seissol::parallel {
template <typename NativeStreamT, typename NativeGraphT, typename NativeInstanceT>
bool SyclNativeGraph<NativeStreamT, NativeGraphT, NativeInstanceT>::canCapture() const {
#if defined(SEISSOL_SYCL_BACKEND_CUDA) || defined(SEISSOL_SYCL_BACKEND_HIP)
  return true;
#else
  return false;
#endif
}
template <typename NativeStreamT, typename NativeGraphT, typename NativeInstanceT>
void SyclNativeGraph<NativeStreamT, NativeGraphT, NativeInstanceT>::beginCapture(StreamT& queue) {
  syclNativeOperation(queue, true, [&](auto& stream) {
    using LocalNativeStreamT = std::decay_t<decltype(stream)>;
    static_assert(std::is_same_v<NativeStreamT, LocalNativeStreamT>, "");
#ifdef SEISSOL_SYCL_BACKEND_CUDA
    if constexpr (std::is_same_v<NativeStreamT, cudaStream_t>) {
      cudaStreamBeginCapture(stream, cudaStreamCaptureModeThreadLocal);
    }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
    if constexpr (std::is_same_v<NativeStreamT, hipStream_t>) {
      hipStreamBeginCapture(stream, hipStreamCaptureModeThreadLocal);
    }
#endif
  });
}
template <typename NativeStreamT, typename NativeGraphT, typename NativeInstanceT>
void SyclNativeGraph<NativeStreamT, NativeGraphT, NativeInstanceT>::endCapture(StreamT& queue) {
  syclNativeOperation(queue, true, [&](auto& stream) {
    using LocalNativeStreamT = std::decay_t<decltype(stream)>;
    static_assert(std::is_same_v<NativeStreamT, LocalNativeStreamT>, "");
#ifdef SEISSOL_SYCL_BACKEND_CUDA
    if constexpr (std::is_same_v<NativeStreamT, cudaStream_t>) {
      cudaStreamEndCapture(stream, &graph);
      instance = std::make_optional<cudaGraphExec_t>();
      cudaGraphInstantiate(&instance.value(), graph, nullptr, nullptr, 0);
    }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
    if constexpr (std::is_same_v<NativeStreamT, hipStream_t>) {
      hipStreamEndCapture(stream, &graph);
      instance = std::make_optional<hipGraphExec_t>();
      hipGraphInstantiate(&instance.value(), graph, nullptr, nullptr, 0);
    }
#endif
  });
}
template <typename NativeStreamT, typename NativeGraphT, typename NativeInstanceT>
void SyclNativeGraph<NativeStreamT, NativeGraphT, NativeInstanceT>::runCapture(StreamT& queue) {
  assert(instance.has_value());
  syclNativeOperation(queue, false, [&](auto& stream) {
    using LocalNativeStreamT = std::decay_t<decltype(stream)>;
    static_assert(std::is_same_v<NativeStreamT, LocalNativeStreamT>, "");
#ifdef SEISSOL_SYCL_BACKEND_CUDA
    if constexpr (std::is_same_v<NativeStreamT, cudaStream_t>) {
      cudaGraphLaunch(instance.value(), stream);
    }
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
    if constexpr (std::is_same_v<NativeStreamT, hipStream_t>) {
      hipGraphLaunch(instance.value(), stream);
    }
#endif
  });
}

#ifdef SYCL_EXT_ONEAPI_GRAPH
void SyclOneapiGraph::canCapture() const { return true; }
void SyclOneapiGraph::beginCapture(StreamT& stream) { graph.begin_recording({stream}); }
void SyclOneapiGraph::endCapture(StreamT& stream) {
  graph.end_recording();
  instance =
      std::move(std::optional<sycl::ext::oneapi::experimental::command_graph<
                    sycl::ext::oneapi::experimental::graph_state::executable>>(graph.finalize()));
}
void SyclOneapiGraph::runCapture(StreamT& stream) {
  assert(instance.has_value());
  stream.submit([&](sycl::handler& handler) { handler.ext_oneapi_graph(instance.value()); });
}
#endif

#ifdef SEISSOL_SYCL_BACKEND_CUDA
template class SyclNativeGraph<cudaStream_t, cudaGraph_t, cudaGraphExec_t>;
#endif
#ifdef SEISSOL_SYCL_BACKEND_HIP
template class SyclNativeGraph<hipStream_t, hipGraph_t, hipGraphInstance_t>;
#endif
} // namespace seissol::parallel

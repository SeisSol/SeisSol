// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/RtcGpuDriver.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Codegen.h"
#include "Expr/Cost.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcGpu.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Datafield/Interpolation.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#if defined(SEISSOL_EXPR_HAVE_NVRTC)
#include <cuda.h>
#include <nvrtc.h>
#endif

#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
#include <hip/hip_runtime.h>
#include <hip/hiprtc.h>
#endif

#if defined(SEISSOL_EXPR_HAVE_SYCL)
#include <sycl/sycl.hpp>
#endif

namespace seissol::expr {

namespace {

namespace df = reader::datafield;
using reader::scripting::DataTable;

/// A loaded kernel, in whichever driver's terms. Held as an opaque pair so the
/// cache and the Kernel below are the same code for both vendors -- the two
/// APIs differ in spelling far more than in shape.
struct DeviceFunction {
  void* module{nullptr};
  void* function{nullptr};
  /// Null when the lowering has no precompute stage.
  void* precompute{nullptr};
  [[nodiscard]] bool valid() const { return function != nullptr; }
};

struct CacheKey {
  std::uint64_t program{0};
  std::uint64_t lowering{0};
  std::uint64_t layout{0};
  GpuTarget target{GpuTarget::Cuda};
  ComputeType type{ComputeType::F64};
  std::string arch;

  bool operator<(const CacheKey& other) const {
    return std::tie(program, lowering, layout, target, type, arch) <
           std::tie(
               other.program, other.lowering, other.layout, other.target, other.type, other.arch);
  }
};

std::mutex& cacheMutex() {
  static std::mutex mutex;
  return mutex;
}

std::map<CacheKey, DeviceFunction>& cache() {
  static std::map<CacheKey, DeviceFunction> functions;
  return functions;
}

// --- launch geometry ---------------------------------------------------------

constexpr unsigned DefaultBlockSize = GpuBlockSize;

/// The grid is capped rather than sized to the point count, because the emitted
/// kernel is a grid-stride loop: a short range still uses every thread it is
/// given, and a long one does not need a block per point. The cap keeps a
/// per-face launch from asking for a grid of one block per node.
constexpr unsigned MaxBlocks = 16384;

[[maybe_unused]] unsigned blocksFor(std::size_t count) {
  if (count == 0) {
    return 0;
  }
  const auto needed = static_cast<unsigned>((count + DefaultBlockSize - 1) / DefaultBlockSize);
  return needed < MaxBlocks ? needed : MaxBlocks;
}

// --- vendor layer ------------------------------------------------------------
//
// Everything below is the only part of the device path that cannot be tested
// without hardware. It is deliberately thin: compile a string, load it, launch
// it. The arithmetic, the argument packing and the rejection rules all sit in
// RtcGpu.cpp, where a host test covers them.

#if defined(SEISSOL_EXPR_HAVE_NVRTC)

bool ensureContext() {
  // The driver API needs a current context, and SeisSol's device layer creates
  // one through the runtime API. Retaining the primary context here rather than
  // creating our own keeps the module in the same context as the buffers the
  // kernel will read.
  static bool initialised = false;
  if (!initialised) {
    if (cuInit(0) != CUDA_SUCCESS) {
      return false;
    }
    initialised = true;
  }
  CUcontext current = nullptr;
  if (cuCtxGetCurrent(&current) == CUDA_SUCCESS && current != nullptr) {
    return true;
  }
  CUdevice device = 0;
  CUcontext primary = nullptr;
  if (cuDeviceGet(&device, 0) != CUDA_SUCCESS ||
      cuDevicePrimaryCtxRetain(&primary, device) != CUDA_SUCCESS ||
      cuCtxSetCurrent(primary) != CUDA_SUCCESS) {
    return false;
  }
  return true;
}

DeviceFunction
    compileCuda(const std::string& source, const std::string& arch, const std::string& dumpName) {
  DeviceFunction loaded;
  if (!ensureContext()) {
    logWarning() << "expr: no usable CUDA context; not using the compiled device backend.";
    return loaded;
  }

  nvrtcProgram program = nullptr;
  if (nvrtcCreateProgram(&program, source.c_str(), "seissol_expr.cu", 0, nullptr, nullptr) !=
      NVRTC_SUCCESS) {
    return loaded;
  }

  // --fmad is left at its default, i.e. contraction ON. Unlike the CPU backend,
  // bitwise agreement with the interpreter is NOT the bar here: turning FMA off
  // costs real throughput on a GPU, and the accepted difference is a documented
  // one-ulp-per-contraction bound rather than none at all.
  const std::string archOption = "--gpu-architecture=compute_" + arch;
  const std::vector<const char*> optionList = {
      archOption.c_str(), "--std=c++17", "-default-device"};

  const nvrtcResult compiled =
      nvrtcCompileProgram(program, static_cast<int>(optionList.size()), optionList.data());
  if (compiled != NVRTC_SUCCESS) {
    std::size_t logSize = 0;
    nvrtcGetProgramLogSize(program, &logSize);
    std::string log(logSize, '\0');
    nvrtcGetProgramLog(program, log.data());
    // The log is the only thing that makes this diagnosable, and a generated
    // source is exactly the case where "it did not compile" is useless on its
    // own.
    logWarning() << "expr: NVRTC rejected the generated kernel:\n" << log.c_str();
    nvrtcDestroyProgram(&program);
    return loaded;
  }

  std::size_t ptxSize = 0;
  nvrtcGetPTXSize(program, &ptxSize);
  std::string ptx(ptxSize, '\0');
  nvrtcGetPTX(program, ptx.data());
  nvrtcDestroyProgram(&program);
  codegen::dump(dumpName + ".ptx", ptx);

  CUmodule module = nullptr;
  if (cuModuleLoadData(&module, ptx.c_str()) != CUDA_SUCCESS) {
    logWarning() << "expr: could not load the generated PTX module.";
    return loaded;
  }
  CUfunction function = nullptr;
  // extern "C" in the emitted source, so the name is not mangled and no
  // nvrtcAddNameExpression round trip is needed.
  if (cuModuleGetFunction(&function, module, "seissol_expr_run") != CUDA_SUCCESS) {
    cuModuleUnload(module);
    logWarning() << "expr: the generated module has no seissol_expr_run.";
    return loaded;
  }
  loaded.module = module;
  loaded.function = function;
  // Optional: only emitted when the lowering has a precompute stage.
  CUfunction precompute = nullptr;
  if (cuModuleGetFunction(&precompute, module, "seissol_expr_precompute") == CUDA_SUCCESS) {
    loaded.precompute = precompute;
  }
  return loaded;
}

void* allocCuda(std::size_t bytes) {
  CUdeviceptr pointer = 0;
  if (bytes == 0 || cuMemAlloc(&pointer, bytes) != CUDA_SUCCESS) {
    return nullptr;
  }
  // Zeroed, which is what makes hoisting work with no further ceremony: the
  // precompute stage writes every hoisted slot before anything reads one, and
  // zero-initialised declared state needs nothing more.
  cuMemsetD8(pointer, 0, bytes);
  return reinterpret_cast<void*>(pointer);
}

void freeCuda(void* pointer) {
  if (pointer != nullptr) {
    cuMemFree(reinterpret_cast<CUdeviceptr>(pointer));
  }
}

bool launchCuda(void* entry, GpuArguments& args, std::size_t count, void* stream) {
  const unsigned blocks = blocksFor(count);
  if (blocks == 0 || entry == nullptr) {
    return blocks == 0;
  }
  const CUresult result = cuLaunchKernel(static_cast<CUfunction>(entry),
                                         blocks,
                                         1,
                                         1,
                                         DefaultBlockSize,
                                         1,
                                         1,
                                         0,
                                         static_cast<CUstream>(stream),
                                         args.data(),
                                         nullptr);
  return result == CUDA_SUCCESS;
}

bool copyToCuda(void* target, const void* source, std::size_t bytes) {
  return ensureContext() &&
         cuMemcpyHtoD(reinterpret_cast<CUdeviceptr>(target), source, bytes) == CUDA_SUCCESS;
}

bool copyFromCuda(void* target, const void* source, std::size_t bytes) {
  return ensureContext() &&
         cuMemcpyDtoH(target, reinterpret_cast<CUdeviceptr>(source), bytes) == CUDA_SUCCESS;
}

/// The work of all streams of the device, which those of SeisSol are, done.
bool synchronizeCuda() { return ensureContext() && cuCtxSynchronize() == CUDA_SUCCESS; }

bool deviceAccessibleCuda(const void* pointer) {
  // Checked rather than declared. A consumer that got this wrong would not get
  // a diagnostic but a fault inside the kernel, and the pointer is the one
  // thing the runtime can answer authoritatively.
  CUmemorytype kind{};
  const CUresult result = cuPointerGetAttribute(
      &kind, CU_POINTER_ATTRIBUTE_MEMORY_TYPE, reinterpret_cast<CUdeviceptr>(pointer));
  if (result != CUDA_SUCCESS) {
    return false;
  }
  return kind == CU_MEMORYTYPE_DEVICE || kind == CU_MEMORYTYPE_UNIFIED;
}

#endif // SEISSOL_EXPR_HAVE_NVRTC

#if defined(SEISSOL_EXPR_HAVE_HIPRTC)

DeviceFunction compileHip(const std::string& source, const std::string& arch) {
  DeviceFunction loaded;

  hiprtcProgram program = nullptr;
  if (hiprtcCreateProgram(&program, source.c_str(), "seissol_expr.hip", 0, nullptr, nullptr) !=
      HIPRTC_SUCCESS) {
    return loaded;
  }

  // Without an architecture, HIPRTC compiles for the current device, which is the one the kernel
  // runs on.
  const std::string archOption = "--offload-arch=" + arch;
  std::vector<const char*> optionList = {"--std=c++17"};
  if (!arch.empty()) {
    optionList.push_back(archOption.c_str());
  }

  const hiprtcResult compiled =
      hiprtcCompileProgram(program, static_cast<int>(optionList.size()), optionList.data());
  if (compiled != HIPRTC_SUCCESS) {
    std::size_t logSize = 0;
    hiprtcGetProgramLogSize(program, &logSize);
    std::string log(logSize, '\0');
    hiprtcGetProgramLog(program, log.data());
    logWarning() << "expr: HIPRTC rejected the generated kernel:\n" << log.c_str();
    hiprtcDestroyProgram(&program);
    return loaded;
  }

  std::size_t codeSize = 0;
  hiprtcGetCodeSize(program, &codeSize);
  std::string code(codeSize, '\0');
  hiprtcGetCode(program, code.data());
  hiprtcDestroyProgram(&program);

  hipModule_t module = nullptr;
  if (hipModuleLoadData(&module, code.data()) != hipSuccess) {
    logWarning() << "expr: could not load the generated HIP module.";
    return loaded;
  }
  hipFunction_t function = nullptr;
  if (hipModuleGetFunction(&function, module, "seissol_expr_run") != hipSuccess) {
    hipModuleUnload(module);
    logWarning() << "expr: the generated module has no seissol_expr_run.";
    return loaded;
  }
  loaded.module = module;
  loaded.function = function;
  hipFunction_t precompute = nullptr;
  if (hipModuleGetFunction(&precompute, module, "seissol_expr_precompute") == hipSuccess) {
    loaded.precompute = precompute;
  }
  return loaded;
}

void* allocHip(std::size_t bytes) {
  void* pointer = nullptr;
  if (bytes == 0 || hipMalloc(&pointer, bytes) != hipSuccess) {
    return nullptr;
  }
  hipMemset(pointer, 0, bytes);
  return pointer;
}

void freeHip(void* pointer) {
  if (pointer != nullptr) {
    hipFree(pointer);
  }
}

bool launchHip(void* entry, GpuArguments& args, std::size_t count, void* stream) {
  const unsigned blocks = blocksFor(count);
  if (blocks == 0 || entry == nullptr) {
    return blocks == 0;
  }
  return hipModuleLaunchKernel(static_cast<hipFunction_t>(entry),
                               blocks,
                               1,
                               1,
                               DefaultBlockSize,
                               1,
                               1,
                               0,
                               static_cast<hipStream_t>(stream),
                               args.data(),
                               nullptr) == hipSuccess;
}

bool copyToHip(void* target, const void* source, std::size_t bytes) {
  return hipMemcpy(target, source, bytes, hipMemcpyHostToDevice) == hipSuccess;
}

bool copyFromHip(void* target, const void* source, std::size_t bytes) {
  return hipMemcpy(target, source, bytes, hipMemcpyDeviceToHost) == hipSuccess;
}

bool synchronizeHip() { return hipDeviceSynchronize() == hipSuccess; }

bool deviceAccessibleHip(const void* pointer) {
  hipPointerAttribute_t attributes{};
  if (hipPointerGetAttributes(&attributes, pointer) != hipSuccess) {
    return false;
  }
  return attributes.type == hipMemoryTypeDevice || attributes.type == hipMemoryTypeUnified;
}

#endif // SEISSOL_EXPR_HAVE_HIPRTC

#if defined(SEISSOL_EXPR_HAVE_SYCL)

namespace syclex = sycl::ext::oneapi::experimental;

/// Build the OpenCL C source into a sycl::kernel on `queue`'s context.
///
/// WHY THE SOURCE IS OpenCL C AND THE ENQUEUE IS SYCL. Raw OpenCL would need
/// sycl::get_native<backend::opencl>() on the queue, which fails when the SYCL
/// backend is Level Zero -- the normal case on an Intel GPU. Going through
/// kernel_compiler keeps the enqueue backend-agnostic. Keeping the SOURCE in
/// OpenCL C rather than SYCL keeps the half that would be expensive to change
/// out of an extension whose own specification says shipping software should
/// not depend on it.
///
/// WHAT IS QUERIED AND WHY. Whether a backend accepts OpenCL C source is a
/// property of the BACKEND, not of the build: a Level Zero device may decline
/// it even where the extension is present. So it is asked at run time and a
/// refusal falls back like any other, rather than being assumed at configure
/// time and failing on the machine that matters.
DeviceFunction compileOpenCl(const std::string& source, void* queuePointer) {
  DeviceFunction loaded;
  if (queuePointer == nullptr) {
    logWarning() << "expr: the OpenCL backend needs a sycl::queue in "
                    "BackendOptions::deviceQueue; using another backend.";
    return loaded;
  }
  auto* queue = static_cast<sycl::queue*>(queuePointer);

  try {
    if (!queue->get_device().ext_oneapi_can_compile(syclex::source_language::opencl)) {
      logWarning() << "expr: this device does not compile OpenCL C at run time; using another "
                      "backend.";
      return loaded;
    }
    auto bundle = syclex::create_kernel_bundle_from_source(
        queue->get_context(), syclex::source_language::opencl, source);
    auto executable = syclex::build(bundle);
    if (!executable.ext_oneapi_has_kernel("seissol_expr_run")) {
      logWarning() << "expr: the generated OpenCL module has no seissol_expr_run.";
      return loaded;
    }
    loaded.module = new sycl::kernel_bundle<sycl::bundle_state::executable>(executable);
    loaded.function = new sycl::kernel(executable.ext_oneapi_get_kernel("seissol_expr_run"));
    if (executable.ext_oneapi_has_kernel("seissol_expr_precompute")) {
      loaded.precompute =
          new sycl::kernel(executable.ext_oneapi_get_kernel("seissol_expr_precompute"));
    }
  } catch (const sycl::exception& error) {
    // The message carries the compiler's own diagnostics, which for generated
    // source is the only thing that makes a failure actionable.
    logWarning() << "expr: the OpenCL kernel did not build:" << error.what();
    return DeviceFunction{};
  }
  return loaded;
}

bool launchOpenCl(void* entry, GpuArguments& args, std::size_t count, void* queuePointer) {
  if (entry == nullptr || queuePointer == nullptr) {
    return false;
  }
  if (count == 0) {
    return true;
  }
  auto* queue = static_cast<sycl::queue*>(queuePointer);
  auto* kernel = static_cast<sycl::kernel*>(entry);

  const std::size_t local = DefaultBlockSize;
  const std::size_t groups = blocksFor(count);
  const std::size_t global = groups * local;

  try {
    queue->submit([&](sycl::handler& handler) {
      // raw_kernel_arg, because the kernel's arity is a property of the PROGRAM
      // and there is no compile-time type list to hand to set_args. This is the
      // flat view of the same packing cuLaunchKernel gets as one struct.
      for (std::size_t i = 0; i < args.fieldCount(); ++i) {
        handler.set_arg(static_cast<int>(i),
                        syclex::raw_kernel_arg(args.fieldData(i), args.fieldSize(i)));
      }
      handler.parallel_for(sycl::nd_range<1>{sycl::range<1>{global}, sycl::range<1>{local}},
                           *kernel);
    });
  } catch (const sycl::exception& error) {
    logError() << "expr: launching the OpenCL kernel failed:" << error.what();
    return false;
  }
  return true;
}

void* allocOpenCl(std::size_t bytes, void* queuePointer) {
  if (bytes == 0 || queuePointer == nullptr) {
    return nullptr;
  }
  auto* queue = static_cast<sycl::queue*>(queuePointer);
  void* pointer = sycl::malloc_device(bytes, *queue);
  if (pointer != nullptr) {
    queue->memset(pointer, 0, bytes).wait();
  }
  return pointer;
}

void freeOpenCl(void* pointer, void* queuePointer) {
  if (pointer != nullptr && queuePointer != nullptr) {
    sycl::free(pointer, *static_cast<sycl::queue*>(queuePointer));
  }
}

bool copyOpenCl(void* target, const void* source, std::size_t bytes, void* queuePointer) {
  if (queuePointer == nullptr) {
    return false;
  }
  try {
    static_cast<sycl::queue*>(queuePointer)->memcpy(target, source, bytes).wait();
  } catch (const sycl::exception& error) {
    logWarning() << "expr: a copy to or from the device failed:" << error.what();
    return false;
  }
  return true;
}

bool synchronizeOpenCl(void* queuePointer) {
  if (queuePointer == nullptr) {
    return false;
  }
  try {
    static_cast<sycl::queue*>(queuePointer)->wait();
  } catch (const sycl::exception& error) {
    logWarning() << "expr: waiting for the device failed:" << error.what();
    return false;
  }
  return true;
}

#endif // SEISSOL_EXPR_HAVE_SYCL

DeviceFunction compileFor(GpuTarget target,
                          const std::string& source,
                          const std::string& arch,
                          void* queue,
                          const std::string& dumpName) {
  codegen::dump(dumpName + (target == GpuTarget::OpenCl ? ".cl" : ".cu"), source);
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return compileCuda(source, arch, dumpName);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return compileHip(source, arch);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    return compileOpenCl(source, queue);
  }
#endif
  static_cast<void>(target);
  static_cast<void>(source);
  static_cast<void>(arch);
  static_cast<void>(queue);
  static_cast<void>(dumpName);
  return {};
}

bool launchFor(GpuTarget target, void* entry, GpuArguments& args, std::size_t count, void* stream) {
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return launchCuda(entry, args, count, stream);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return launchHip(entry, args, count, stream);
  }
#endif
  static_cast<void>(target);
  static_cast<void>(entry);
  static_cast<void>(args);
  static_cast<void>(count);
  static_cast<void>(stream);
  return false;
}

void* allocFor(GpuTarget target, std::size_t bytes, void* queue) {
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    return allocOpenCl(bytes, queue);
  }
#endif
  static_cast<void>(queue);
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return allocCuda(bytes);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return allocHip(bytes);
  }
#endif
  static_cast<void>(target);
  static_cast<void>(bytes);
  return nullptr;
}

void freeFor(GpuTarget target, void* pointer, void* queue) {
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    freeOpenCl(pointer, queue);
    return;
  }
#endif
  static_cast<void>(queue);
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    freeCuda(pointer);
    return;
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    freeHip(pointer);
    return;
  }
#endif
  static_cast<void>(target);
  static_cast<void>(pointer);
}

/// A copy of `bytes` from the host to the device; done when it returns.
bool copyToDeviceFor(
    GpuTarget target, void* device, const void* host, std::size_t bytes, void* queue) {
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    return copyOpenCl(device, host, bytes, queue);
  }
#endif
  static_cast<void>(queue);
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return copyToCuda(device, host, bytes);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return copyToHip(device, host, bytes);
  }
#endif
  static_cast<void>(target);
  static_cast<void>(device);
  static_cast<void>(host);
  static_cast<void>(bytes);
  return false;
}

/// A copy of `bytes` from the device to the host; done when it returns.
bool copyFromDeviceFor(
    GpuTarget target, void* host, const void* device, std::size_t bytes, void* queue) {
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    return copyOpenCl(host, device, bytes, queue);
  }
#endif
  static_cast<void>(queue);
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return copyFromCuda(host, device, bytes);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return copyFromHip(host, device, bytes);
  }
#endif
  static_cast<void>(target);
  static_cast<void>(device);
  static_cast<void>(host);
  static_cast<void>(bytes);
  return false;
}

/// The work done that was submitted before: on the whole device for CUDA and HIP, on `queues` for
/// OpenCL, whose kernels each run on the queue they were built for.
bool synchronizeFor(GpuTarget target, const std::vector<void*>& queues) {
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    bool done = true;
    for (void* queue : queues) {
      done = synchronizeOpenCl(queue) && done;
    }
    return done;
  }
#endif
  static_cast<void>(queues);
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return synchronizeCuda();
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return synchronizeHip();
  }
#endif
  static_cast<void>(target);
  return false;
}

bool (*accessibilityCheck(GpuTarget target))(const void*) {
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return deviceAccessibleCuda;
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return deviceAccessibleHip;
  }
#endif
  static_cast<void>(target);
  return nullptr;
}

// --- grids on the device -----------------------------------------------------

/// The copy of a grid on the device that the kernels sample (a GridStore::Mirror, refreshed where
/// the grid changes): its array, and its descriptor (codegen::GridDescriptor) at an address that
/// stays the same, through which a kernel reads the array -- also one recorded in a graph, which
/// then samples the window of a time-dependent grid that is resident when it runs.
///
/// A refresh copies the whole array, after the work on the device is done: the store refreshes
/// at a synchronization point only, where the window of a time-dependent grid has moved.
class DeviceGridMirror final : public df::GridStore::Mirror {
  public:
  /// The kind of copy in the store (GridStore::mirror()), one per target.
  static int keyOf(GpuTarget target) { return 0x45585000 + static_cast<int>(target); }

  DeviceGridMirror(GpuTarget target, void* queue, const df::GridDesc& desc)
      : target_(target), queue_(queue), path_(desc.path),
        fails_(desc.boundary == df::BoundaryMode::Fail) {
    addQueue(queue);
  }

  // Frees nothing: the store that owns it may live until the process ends, past the runtime of the
  // device. finish() frees the memory, at the end of the run.
  ~DeviceGridMirror() override = default;

  DeviceGridMirror(const DeviceGridMirror&) = delete;
  DeviceGridMirror& operator=(const DeviceGridMirror&) = delete;
  DeviceGridMirror(DeviceGridMirror&&) = delete;
  DeviceGridMirror& operator=(DeviceGridMirror&&) = delete;

  /// The descriptor on the device, or null where the copy failed.
  [[nodiscard]] void* descriptor() const { return valid_ ? descriptor_ : nullptr; }
  [[nodiscard]] std::size_t bytes() const { return bytes_; }

  /// A queue a kernel that samples the copy runs on, which a refresh waits for (OpenCL).
  void addQueue(void* queue) {
    if (std::find(queues_.begin(), queues_.end(), queue) == queues_.end()) {
      queues_.push_back(queue);
    }
  }

  void refresh(const df::Grid& grid) override {
    if (finished_ || (descriptor_ != nullptr && revision_ == grid.revision())) {
      return;
    }
    const df::ArrayView* view = grid.arrayView();
    if (view == nullptr) {
      logError() << "expr: the grid" << path_ << "on the device no longer keeps an array.";
      return;
    }
    if (descriptor_ != nullptr) {
      // the kernels that read the copy are done with it; and what they found is reported
      synchronizeFor(target_, queues_);
      report();
    }
    valid_ = false;
    const std::size_t bytes = codegen::gridBytes(*view);
    if (bytes != bytes_) {
      freeFor(target_, data_, queue_);
      data_ = allocFor(target_, bytes, queue_);
      bytes_ = data_ != nullptr ? bytes : 0;
    }
    if (bytes > 0 &&
        (data_ == nullptr || !copyToDeviceFor(target_, data_, view->data, bytes, queue_))) {
      logWarning() << "expr: could not copy the" << bytes << "bytes of the grid" << path_
                   << "to the device.";
      return;
    }
    if (descriptor_ == nullptr) {
      descriptor_ = allocFor(target_, sizeof(codegen::GridDescriptor), queue_);
    }
    const codegen::GridDescriptor image = codegen::gridDescriptor(*view, data_);
    if (descriptor_ == nullptr ||
        !copyToDeviceFor(target_, descriptor_, &image, sizeof(image), queue_)) {
      logWarning() << "expr: could not copy the descriptor of the grid" << path_
                   << "to the device.";
      return;
    }
    // done for a kernel on any stream that follows
    synchronizeFor(target_, queues_);
    revision_ = grid.revision();
    valid_ = true;
  }

  void finish(const df::Grid& /*grid*/) override {
    if (finished_) {
      return;
    }
    finished_ = true;
    if (descriptor_ == nullptr) {
      return;
    }
    synchronizeFor(target_, queues_);
    report();
    freeFor(target_, data_, queue_);
    freeFor(target_, descriptor_, queue_);
    data_ = nullptr;
    descriptor_ = nullptr;
    valid_ = false;
  }

  private:
  /// A lookup outside the volume, of a grid whose boundary mode fails there: an error, as on the
  /// host, which a kernel cannot raise itself (codegen::GridDescriptor::outside).
  void report() const {
    if (!fails_ || !valid_) {
      return;
    }
    std::uint64_t outside = 0;
    const auto* field =
        static_cast<const char*>(descriptor_) + offsetof(codegen::GridDescriptor, outside);
    if (copyFromDeviceFor(target_, &outside, field, sizeof(outside), queue_) && outside != 0) {
      logError() << "A kernel on the device sampled the grid" << path_
                 << "outside its volume, where its boundary mode is to fail.";
    }
  }

  GpuTarget target_;
  void* queue_;
  std::vector<void*> queues_;
  std::string path_;
  bool fails_;
  void* data_{nullptr};
  std::size_t bytes_{0};
  void* descriptor_{nullptr};
  std::uint64_t revision_{0};
  bool valid_{false};
  bool finished_{false};
};

/// The descriptors on the device of the grids `program` looks up, by GridId (null for one it does
/// not look up), with the types of their elements in `layout`; copies a grid to the device where
/// none is there yet. False where a grid cannot go there, which warns.
bool deviceGrids(const Program& program,
                 const LoweredProgram& lowered,
                 GpuTarget target,
                 void* queue,
                 df::GridStore& grids,
                 const std::vector<std::size_t>& gridIds,
                 GpuLayout& layout,
                 std::vector<void*>& pointers) {
  const std::size_t count = program.grids().size();
  // the coordinates and the components each grid is sampled with
  std::vector<std::int32_t> dimensions(count, 0);
  std::vector<std::int32_t> components(count, 0);
  for (std::size_t grid = 0; grid < count && grid < lowered.gridUses().size(); ++grid) {
    dimensions[grid] = lowered.gridUses()[grid].dimensions;
    components[grid] = lowered.gridUses()[grid].components;
  }
  layout.grids.assign(count, df::ElementType::Double);
  pointers.assign(count, nullptr);
  for (std::size_t grid = 0; grid < count; ++grid) {
    if (dimensions[grid] == 0) {
      continue;
    }
    const auto& desc = program.grids()[grid];
    const auto reject = [&](const char* why) {
      logWarning() << "expr: this program cannot run on a device because"
                   << describe(GpuRejection::Lookup) << "-- the grid" << desc.path << why
                   << "-- using another backend.";
      return false;
    };
    if (grid >= gridIds.size()) {
      return reject("is not in the store");
    }
    const std::size_t id = gridIds[grid];
    const df::Grid& source = grids.get(id);
    const df::ArrayView* view = source.arrayView();
    if (view == nullptr) {
      return reject("keeps no array");
    }
    if (static_cast<std::int32_t>(view->dimensions) != dimensions[grid] ||
        components[grid] > static_cast<std::int32_t>(view->components) ||
        view->dimensions > df::MaxGridDimensions) {
      return reject("has other dimensions or fewer components than the lookups");
    }
    layout.grids[grid] = view->type;
    auto* mirror =
        dynamic_cast<DeviceGridMirror*>(grids.mirror(id, DeviceGridMirror::keyOf(target)));
    if (mirror == nullptr) {
      mirror = &dynamic_cast<DeviceGridMirror&>(
          grids.attachMirror(id,
                             DeviceGridMirror::keyOf(target),
                             std::make_unique<DeviceGridMirror>(target, queue, desc)));
      mirror->refresh(source);
      if (mirror->descriptor() != nullptr) {
        logInfo() << "expr: the grid" << desc.path << "is copied to the device (" << mirror->bytes()
                  << "bytes).";
      }
    }
    mirror->addQueue(queue);
    if (mirror->descriptor() == nullptr) {
      return reject("could not be copied to the device");
    }
    pointers[grid] = mirror->descriptor();
  }
  return true;
}

// --- kernel ------------------------------------------------------------------

class RtcGpuKernel final : public Kernel {
  public:
  RtcGpuKernel(Binding& binding,
               LoweredProgram lowered,
               DeviceFunction function,
               GpuTarget target,
               std::size_t elementWidth,
               void* queue,
               std::vector<void*> grids)
      : binding_(&binding), lowered_(std::move(lowered)), function_(function), target_(target),
        queue_(queue), grids_(std::move(grids)) {
    // The hoisted values, and the states the table does not keep; both zeroed by allocFor, which
    // is all a zero-initialised state needs (gpuRejection refuses the others).
    const auto stateSlots = static_cast<std::size_t>(lowered_.stateSlotCount());
    const auto hoisted = static_cast<std::size_t>(lowered_.persistentSlotCount()) - stateSlots;
    persistent_ = allocate(hoisted * binding.numPoints() * elementWidth, "hoisted values");
    bool ownsState = false;
    for (std::size_t i = 0; i < binding.states().size(); ++i) {
      ownsState |= !binding.stateKept(i);
    }
    if (ownsState) {
      ownState_ = allocate(stateSlots * binding.numPoints() * elementWidth, "states");
    }
  }

  ~RtcGpuKernel() override {
    freeFor(target_, persistent_, queue_);
    freeFor(target_, ownState_, queue_);
  }

  RtcGpuKernel(const RtcGpuKernel&) = delete;
  RtcGpuKernel& operator=(const RtcGpuKernel&) = delete;
  RtcGpuKernel(RtcGpuKernel&&) = delete;
  RtcGpuKernel& operator=(RtcGpuKernel&&) = delete;

  /// Fill the hoisted slots. Safe to call again, and MEANT to be: anything that
  /// changes an input LowerOptions declared invariant -- a velocity model
  /// swapped mid-run by the instantaneous time mirroring, say -- makes the
  /// hoisted values stale, and re-running this is the whole remedy. The buffer
  /// is not reallocated, so the call costs one launch and nothing else.
  void precompute(const DataTable& /*table*/) override {
    if (function_.precompute == nullptr) {
      return;
    }
    KernelArgs args{};
    args.first = 0;
    args.count = binding_->numPoints();
    GpuArguments packed(*binding_, args, persistent_, ownState_, grids_);
    if (!launchFor(target_, function_.precompute, packed, args.count, streamFor(args))) {
      logError() << "expr: launching the device precompute stage failed.";
    }
    precomputed_ = true;
  }

  void run(const DataTable& /*table*/) override {
    // The bound bases are used as they are, which is the whole table. There is
    // no table-specific path on a device: the accessors it would need are
    // std::functions.
    KernelArgs args{};
    args.first = 0;
    args.count = binding_->numPoints();
    run(args);
  }

  void run(const KernelArgs& args) override {
    if (function_.precompute != nullptr && !precomputed_) {
      // Same guard as everywhere else: the hoisted slots would read as the
      // zeros the allocation left, which is a plausible wrong answer rather
      // than a fault.
      logError() << "expr: the kernel has a precompute stage that was never run; call "
                    "Kernel::precompute() from prepare().";
    }
    GpuArguments packed(*binding_, args, persistent_, ownState_, grids_);
    if (!launchFor(target_, function_.function, packed, args.count, streamFor(args))) {
      logError() << "expr: launching the compiled device kernel failed.";
    }
  }

  [[nodiscard]] BackendKind kind() const override {
    switch (target_) {
    case GpuTarget::Cuda:
      return BackendKind::RtcCuda;
    case GpuTarget::Hip:
      return BackendKind::RtcHip;
    case GpuTarget::OpenCl:
      return BackendKind::RtcOpenCl;
    }
    return BackendKind::RtcCuda;
  }

  private:
  [[nodiscard]] void* allocate(std::size_t bytes, const char* what) const {
    if (bytes == 0) {
      return nullptr;
    }
    void* pointer = allocFor(target_, bytes, queue_);
    if (pointer == nullptr) {
      logError() << "expr: could not allocate" << bytes << "bytes of device memory for the" << what
                 << ".";
    }
    return pointer;
  }

  /// The OpenCL path enqueues on the queue it was BUILT against -- the kernel
  /// bundle belongs to that queue's context -- so a per-call stream would have
  /// to come from the same context anyway. Everywhere else the per-call stream
  /// wins, which is what lets a caller overlap two evaluations.
  [[nodiscard]] void* streamFor(const KernelArgs& args) const {
    if (target_ == GpuTarget::OpenCl) {
      return queue_;
    }
    return args.stream;
  }

  Binding* binding_;
  LoweredProgram lowered_;
  DeviceFunction function_;
  GpuTarget target_;
  void* queue_{nullptr};
  /// the descriptors of the grids it samples, on the device (deviceGrids())
  std::vector<void*> grids_;
  void* persistent_{nullptr};
  void* ownState_{nullptr};
  bool precomputed_{false};
};

} // namespace

bool gpuRuntimeAvailable(GpuTarget target) {
#if defined(SEISSOL_EXPR_HAVE_SYCL)
  if (target == GpuTarget::OpenCl) {
    return true;
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return true;
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return true;
  }
#endif
  static_cast<void>(target);
  return false;
}

bool deviceRuntimeAvailable() {
  switch (deviceBackendOfBuild()) {
  case BackendKind::RtcCuda:
    return gpuRuntimeAvailable(GpuTarget::Cuda);
  case BackendKind::RtcHip:
    return gpuRuntimeAvailable(GpuTarget::Hip);
  case BackendKind::RtcOpenCl:
    return gpuRuntimeAvailable(GpuTarget::OpenCl);
  default:
    return false;
  }
}

std::size_t rtcGpuCacheSize() {
  const std::lock_guard<std::mutex> lock(cacheMutex());
  return cache().size();
}

std::unique_ptr<Kernel> makeRtcGpuKernel(const Program& program,
                                         Binding& binding,
                                         LoweredProgram lowered,
                                         const BackendOptions& options,
                                         GpuTarget target,
                                         df::GridStore& grids,
                                         const std::vector<std::size_t>& gridIds) {
  if (!gpuRuntimeAvailable(target)) {
    logWarning() << "expr: this build has no runtime compiler for the requested device backend.";
    return nullptr;
  }

  const GpuRejection rejection =
      gpuRejection(program, lowered, binding, accessibilityCheck(target));
  if (rejection != GpuRejection::None) {
    logWarning() << "expr: this program cannot run on a device because" << describe(rejection)
                 << "-- using another backend.";
    return nullptr;
  }

  GpuLayout layout = gpuLayoutOf(binding);
  // the grids it samples, on the device
  std::vector<void*> gridPointers;
  if (lowered.usesGrids() &&
      !deviceGrids(
          program, lowered, target, options.deviceQueue, grids, gridIds, layout, gridPointers)) {
    return nullptr;
  }
  // PTX for compute_75 runs on every later device, by JIT; a HIP code object has to be for the
  // device itself, which HIPRTC compiles for when it is given none.
  std::string arch = options.arch;
  if (arch.empty() && target == GpuTarget::Cuda) {
    arch = "75";
  }
  const CacheKey key{program.fingerprint(),
                     options.lowering.fingerprint(),
                     layout.fingerprint(),
                     target,
                     program.computeType(),
                     arch};

  DeviceFunction function;
  {
    // Held across the compile. This runs from prepare(), where serialising a
    // handful of compiles costs nothing and two threads building the same
    // kernel would be pure waste.
    const std::lock_guard<std::mutex> lock(cacheMutex());
    auto found = cache().find(key);
    if (found == cache().end()) {
      const std::string source = emitGpuSource(program, lowered, layout, target);
      std::ostringstream dumpName;
      dumpName << "expr-" << std::hex << program.fingerprint() << std::dec << "-"
               << (program.computeType() == ComputeType::F32 ? "f32" : "f64") << "-" << arch;
      found =
          cache()
              .emplace(key, compileFor(target, source, arch, options.deviceQueue, dumpName.str()))
              .first;
    }
    function = found->second;
  }
  if (!function.valid()) {
    return nullptr;
  }

  if (!options.quiet) {
    logInfo()
        << "expr: compiled device kernel --" << lowered.summary().c_str() << "--"
        << cost(program, lowered, program.computeType()).summary(program.computeType()).c_str()
        << "-- for" << arch;
  }
  // the hoisted values and the states the kernel keeps, in the type the data moves in
  const std::size_t width = lowered.dataType() == ComputeType::F32 ? 4 : 8;
  return std::make_unique<RtcGpuKernel>(binding,
                                        std::move(lowered),
                                        function,
                                        target,
                                        width,
                                        options.deviceQueue,
                                        std::move(gridPointers));
}

} // namespace seissol::expr

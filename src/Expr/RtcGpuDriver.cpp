// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/RtcGpuDriver.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Cost.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Expr/RtcGpu.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
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

namespace seissol::expr {

namespace {

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

constexpr unsigned DefaultBlockSize = 256;

/// The grid is capped rather than sized to the point count, because the emitted
/// kernel is a grid-stride loop: a short range still uses every thread it is
/// given, and a long one does not need a block per point. The cap keeps a
/// per-face launch from asking for a grid of one block per node.
constexpr unsigned MaxBlocks = 4096;

unsigned blocksFor(std::size_t count) {
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

DeviceFunction compileCuda(const std::string& source, const std::string& arch) {
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

  const std::string archOption = "--offload-arch=" + arch;
  const std::vector<const char*> optionList = {archOption.c_str(), "--std=c++17"};

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

bool deviceAccessibleHip(const void* pointer) {
  hipPointerAttribute_t attributes{};
  if (hipPointerGetAttributes(&attributes, pointer) != hipSuccess) {
    return false;
  }
  return attributes.type == hipMemoryTypeDevice || attributes.type == hipMemoryTypeUnified;
}

#endif // SEISSOL_EXPR_HAVE_HIPRTC

DeviceFunction compileFor(GpuTarget target, const std::string& source, const std::string& arch) {
#if defined(SEISSOL_EXPR_HAVE_NVRTC)
  if (target == GpuTarget::Cuda) {
    return compileCuda(source, arch);
  }
#endif
#if defined(SEISSOL_EXPR_HAVE_HIPRTC)
  if (target == GpuTarget::Hip) {
    return compileHip(source, arch);
  }
#endif
  static_cast<void>(target);
  static_cast<void>(source);
  static_cast<void>(arch);
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

void* allocFor(GpuTarget target, std::size_t bytes) {
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

void freeFor(GpuTarget target, void* pointer) {
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

// --- kernel ------------------------------------------------------------------

class RtcGpuKernel final : public Kernel {
  public:
  RtcGpuKernel(Binding& binding,
               LoweredProgram lowered,
               DeviceFunction function,
               GpuTarget target,
               std::size_t elementWidth)
      : binding_(&binding), lowered_(std::move(lowered)), function_(function), target_(target) {
    const auto slots = static_cast<std::size_t>(lowered_.persistentSlotCount());
    persistentBytes_ = slots * binding.numPoints() * elementWidth;
    if (persistentBytes_ > 0) {
      persistent_ = allocFor(target_, persistentBytes_);
      if (persistent_ == nullptr) {
        logError() << "expr: could not allocate" << persistentBytes_
                   << "bytes of device memory for the persistent buffer.";
      }
    }
  }

  ~RtcGpuKernel() override { freeFor(target_, persistent_); }

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
    GpuArguments packed(*binding_, args, persistent_);
    if (!launchFor(target_, function_.precompute, packed, args.count, args.stream)) {
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
    GpuArguments packed(*binding_, args, persistent_);
    if (!launchFor(target_, function_.function, packed, args.count, args.stream)) {
      logError() << "expr: launching the compiled device kernel failed.";
    }
  }

  [[nodiscard]] BackendKind kind() const override {
    return target_ == GpuTarget::Cuda ? BackendKind::RtcCuda : BackendKind::RtcHip;
  }

  private:
  Binding* binding_;
  LoweredProgram lowered_;
  DeviceFunction function_;
  GpuTarget target_;
  void* persistent_{nullptr};
  std::size_t persistentBytes_{0};
  bool precomputed_{false};
};

} // namespace

bool gpuRuntimeAvailable(GpuTarget target) {
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

std::size_t rtcGpuCacheSize() {
  const std::lock_guard<std::mutex> lock(cacheMutex());
  return cache().size();
}

std::unique_ptr<Kernel> makeRtcGpuKernel(const Program& program,
                                         Binding& binding,
                                         LoweredProgram lowered,
                                         const BackendOptions& options,
                                         GpuTarget target) {
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

  const GpuLayout layout = gpuLayoutOf(binding);
  const std::string arch = options.arch.empty() ? std::string("70") : options.arch;
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
      found = cache().emplace(key, compileFor(target, source, arch)).first;
    }
    function = found->second;
  }
  if (!function.valid()) {
    return nullptr;
  }

  logInfo() << "expr: compiled device kernel --" << lowered.summary().c_str() << "--"
            << cost(lowered, program.computeType()).summary(program.computeType()).c_str()
            << "-- for" << arch;
  const std::size_t width = program.computeType() == ComputeType::F32 ? 4 : 8;
  return std::make_unique<RtcGpuKernel>(binding, std::move(lowered), function, target, width);
}

} // namespace seissol::expr

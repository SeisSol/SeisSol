// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_RTCGPU_H_
#define SEISSOL_SRC_EXPR_RTCGPU_H_

// Source generation for the device backends.
//
// WHAT IS DIFFERENT FROM THE CPU BACKEND, and it is not the arithmetic. That
// comes from codegen::emitStageBody, the same text RtcCpu emits, from the same
// table the interpreter evaluates.
//
// What differs is the DATA PATH. RtcCpu gathers a tile on the host and hands
// the kernel a dense buffer; a device kernel cannot, because the gather runs
// through DataEntry::accessor, which is a std::function -- not slow on a GPU,
// but impossible. So the emitted kernel reads DataTable::StridedView directly:
// each column arrives as a base pointer plus a byte stride and offset, and the
// element type is baked in at generation time so no per-point type dispatch
// survives into the kernel.
//
// WHY TYPES ARE BAKED BUT STRIDES ARE NOT. A type switch per point would cost
// more than the arithmetic, so it has to disappear at generation time -- which
// means the column TYPES belong in the kernel cache key. Strides and offsets do
// not: they are scalars the kernel reads once per thread, and baking them would
// give one kernel per binding layout instead of one per program, for no
// measurable gain.
//
// WHAT IS DELIBERATELY NOT HERE. The NVRTC/HIPRTC driver layer. This header
// produces a string and answers whether a program can be compiled at all; it
// pulls in no vendor headers, so it builds and is testable everywhere -- and
// the emitted source can be compiled and RUN on the host, which is how its
// arithmetic is checked against the interpreter without a device.

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace seissol::expr {

enum class GpuTarget : std::uint8_t {
  Cuda,
  Hip,
  /// OpenCL C, for Intel.
  ///
  /// The SOURCE language, not the runtime. It is enqueued through SYCL's
  /// kernel_compiler extension on SeisSol's existing sycl::queue rather than
  /// through raw OpenCL, because raw OpenCL would need get_native<opencl>() on
  /// the queue -- which fails on a Level Zero backend, i.e. the normal case on
  /// an Intel GPU.
  ///
  /// OpenCL C rather than SYCL source, although the extension accepts both,
  /// because the extension itself is experimental and its own specification
  /// says shipping software should not depend on it. Keeping the source in the
  /// stable half means that if the extension changes, or is replaced by raw
  /// OpenCL or a SPIR-V path, only the enqueue moves.
  OpenCl
};

/// Whether `target` takes its arguments as one by-value struct (CUDA, HIP) or
/// as a flat parameter list (OpenCL).
///
/// Not a style choice. OpenCL C does not portably allow pointers inside a
/// struct passed as a KERNEL argument, so that target has to spell the
/// parameters out -- and its enqueue is clSetKernelArg per index anyway, which
/// wants exactly that. GpuArguments produces both views from one packing.
[[nodiscard]] bool usesArgumentStruct(GpuTarget target);

/// Why a program cannot go to a device backend. Reported rather than returned
/// as a bare bool, because "fell back to the interpreter" without a reason is
/// the kind of log line nobody acts on.
enum class GpuRejection : std::uint8_t {
  None,
  /// Samples a data grid. A lookup is a batch call into GridSampler, not a
  /// lane-local expression; the device sampler is separate work.
  Lookup,
  /// A column was bound with bindComputed, so there is no address arithmetic
  /// behind it and the kernel has nothing to read.
  ComputedColumn,
  /// The point set had to be reordered to make the groups contiguous, and the
  /// permutation lives in host memory. Refused rather than silently ignored:
  /// evaluating the unpermuted order would put results on the wrong points.
  Permuted,
  /// Declares state with a non-zero initial value.
  ///
  /// NARROWED (reported, Package 5). This used to refuse any persistent slot at
  /// all, which also refused HOISTING -- and hoisting is precisely what a
  /// material-dependent boundary condition needs, since wave speeds are
  /// time-independent and recomputing them every step is most of the work
  /// (measured: 114 weighted ops per point against 54 when hoisted).
  ///
  /// The device buffer is zeroed on allocation and the precompute stage writes
  /// the hoisted slots, so hoisting needs nothing else. What is still missing
  /// is a fill for declared state whose initial value is not zero; that wants a
  /// small kernel and no program in sight uses it.
  StatefulProgram,
  /// At least one base pointer is not device-accessible.
  HostPointer,
  /// A column's stride or offset is not a whole number of its own elements.
  ///
  /// The device kernel loads through a typed pointer, which is undefined on a
  /// misaligned address -- and on a GPU that is a fault or a silently wrong
  /// value, not a slow path. Every DataTable bind form produces element
  /// multiples (bindView by construction, bindMemberView because a member is
  /// aligned within its struct), so this cannot trigger for anything the ABI
  /// itself builds; it is here because the check is two comparisons and the
  /// alternative is an invariant that only holds by argument.
  Misaligned
};

[[nodiscard]] const char* describe(GpuRejection rejection);

/// Whether this program and binding can be evaluated on a device, and why not.
///
/// `deviceAccessible` answers whether a base pointer is reachable from the
/// device. It is a callback rather than a flag on the binding because only the
/// runtime can answer it -- cudaPointerGetAttributes and its HIP equivalent --
/// and asking the consumer to declare it would turn a wrong declaration into a
/// fault inside the kernel. Pass nullptr to skip the check, which is what a
/// codegen test wants.
[[nodiscard]] GpuRejection gpuRejection(const Program& program,
                                        const LoweredProgram& lowered,
                                        const Binding& binding,
                                        bool (*deviceAccessible)(const void*));

/// One column's element type and role, in Binding slot order. Part of the
/// kernel cache key: the types are baked into the emitted source.
struct GpuLayout {
  std::vector<reader::scripting::DataType> inputs;
  std::vector<reader::scripting::DataType> outputs;

  [[nodiscard]] std::uint64_t fingerprint() const;
};

[[nodiscard]] GpuLayout gpuLayoutOf(const Binding& binding);

/// Kernel parameters, packed in the order emitGpuSource declares them.
///
/// This exists so the packing order is stated ONCE. The emitted signature is
/// generated, so its arity depends on the program -- a driver that spelled the
/// order out by hand would be a second description of a generated thing, and
/// the failure mode is silent: parameters land one slot over, the count becomes
/// a pointer, and the kernel writes nothing at all rather than crashing.
///
/// The layout matches cuLaunchKernel's kernelParams: each entry POINTS AT an
/// argument value, and the values live in this object, so it has to outlive the
/// launch. Since the kernel takes a single by-value struct, there is exactly
/// one entry and this object holds that struct's byte image.
class GpuArguments {
  public:
  /// `persistent` may be null when the lowering has no persistent slots, which
  /// gpuRejection() already guarantees for a program that reaches a device.
  GpuArguments(const Binding& binding, const KernelArgs& args, void* persistent);

  /// The cuLaunchKernel view: one entry, pointing at the argument struct.
  [[nodiscard]] void** data() { return pointers_.data(); }
  [[nodiscard]] std::size_t size() const { return pointers_.size(); }

  /// The clSetKernelArg view: one entry per parameter, in declaration order.
  /// Same bytes, same order -- the two views cannot disagree because there is
  /// one packing behind them.
  [[nodiscard]] std::size_t fieldCount() const { return fields_.size(); }
  [[nodiscard]] const void* fieldData(std::size_t index) const {
    return image_.data() + fields_[index].offset;
  }
  [[nodiscard]] std::size_t fieldSize(std::size_t index) const { return fields_[index].size; }

  private:
  void append(const void* value, std::size_t bytes);

  /// The by-value struct's byte image. pointers_ points into it, so it is
  /// reserved up front and never grown afterwards.
  struct Field {
    std::size_t offset{0};
    std::size_t size{0};
  };

  std::vector<unsigned char> image_;
  std::vector<Field> fields_;
  std::vector<void*> pointers_;
};

/// The kernel source. Emits `seissol_expr_run`, and `seissol_expr_precompute`
/// only when the lowering has a precompute stage.
[[nodiscard]] std::string emitGpuSource(const Program& program,
                                        const LoweredProgram& lowered,
                                        const GpuLayout& layout,
                                        GpuTarget target);

/// A `void**`-taking wrapper around the kernel, generated from the same layout
/// as the signature.
///
/// The driver launches through cuLaunchKernel, which unpacks the array itself;
/// this is for anything that has to CALL the kernel as an ordinary function --
/// notably the test that compiles the emitted source for the host and checks it
/// against the interpreter. Generating it from the layout is the point: a
/// hand-written caller would encode the argument order a second time, which is
/// exactly what GpuArguments exists to avoid.
[[nodiscard]] std::string emitGpuHostTrampoline(const GpuLayout& layout,
                                                const std::string& computeType);

/// The same, for the flat-parameter targets: takes an array of POINTERS to the
/// parameter values, which is what GpuArguments' field view hands out and what
/// clSetKernelArg consumes one at a time.
[[nodiscard]] std::string emitGpuHostTrampolineFlat(const GpuLayout& layout,
                                                    const std::string& computeType);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_RTCGPU_H_

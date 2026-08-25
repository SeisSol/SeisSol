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

enum class GpuTarget : std::uint8_t { Cuda, Hip };

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
  /// Carries state or hoisted values. The persistent buffer is a host
  /// allocation owned by Binding; giving it a device counterpart is separate
  /// work, and reading host memory from the kernel would fault.
  PersistentState,
  /// At least one base pointer is not device-accessible.
  HostPointer
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
[[nodiscard]] GpuRejection gpuRejection(const LoweredProgram& lowered,
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

  [[nodiscard]] void** data() { return pointers_.data(); }
  [[nodiscard]] std::size_t size() const { return pointers_.size(); }

  private:
  void append(const void* value, std::size_t bytes);

  /// The by-value struct's byte image. pointers_ points into it, so it is
  /// reserved up front and never grown afterwards.
  std::vector<unsigned char> image_;
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

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_RTCGPU_H_

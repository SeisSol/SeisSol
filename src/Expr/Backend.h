// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_BACKEND_H_
#define SEISSOL_SRC_EXPR_BACKEND_H_

// Backends turn a (Program, Binding) pair into something callable. Everything
// above this line is backend-agnostic; everything below it is free to be as
// specialised as it likes, because the tile layout is fixed by Binding.
//
// Note on reuse from the derived-output side: rtc.hpp and cache.hpp are reusable
// as MECHANISM but not as-is. CpuProgram hardcodes three entry points
// (init_state / perstep_cell / finalize_cell) that only make sense for the
// stateful temporal-reduction kernel; the scripting path wants a single
// stateless `evaluate`. Both need to be parameterised on an entry-point list.
// KernelCache is generic in `Program` but its key derivation goes through
// ir_hash(Arena, Plan) — that becomes Program::fingerprint() here.

#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <memory>
#include <string>

namespace seissol::expr {

enum class BackendKind : std::uint8_t {
  Interpreter, // tiled batch/SoA walk; always available, the correctness oracle
  RtcCpu,      // host compiler -> shared object
  RtcCuda,     // NVRTC
  RtcHip,      // hiprtc
  Texture,     // NOT IMPLEMENTED — hardware-sampled grids, see stub below
  Distributed  // NOT IMPLEMENTED — sharded grids, see stub below
};

const char* name(BackendKind kind);

/// One call's worth of raw arguments, in Binding slot order.
struct KernelArgs {
  /// inputs[i] is the base for Program::inputs()[i], outputs[j] for outputs()[j].
  /// A null entry means "use the base the Binding was bound with", so a caller
  /// that only moves some of the columns need not restate the others.
  const void* const* inputs{nullptr};
  void* const* outputs{nullptr};
  std::size_t inputCount{0};
  std::size_t outputCount{0};

  /// Point range, in the Binding's (possibly permuted) index space.
  /// count == 1 is the element-wise form; nothing about the kernel changes.
  std::size_t first{0};
  std::size_t count{0};

  /// Opaque queue/stream, matching SeisSol's existing void* convention
  /// (Layer::synchronizeTo and friends). Null means synchronous on the host.
  /// Ignored by the interpreter and by RtcCpu.
  void* stream{nullptr};

  /// Whether the backend may open a parallel region.
  ///
  /// Default OFF, and that is a correctness default rather than a conservative
  /// one: the element-wise call site is a loop over faces that is itself
  /// already inside an OpenMP region, and nesting one there would be wrong
  /// regardless of how it performs. The dense run(table) path sets it.
  bool allowThreads{false};
};

// A prepared, ready-to-run instance of one Program against one Binding.
class Kernel {
  public:
  Kernel() = default;
  virtual ~Kernel() = default;
  Kernel(const Kernel&) = delete;
  Kernel& operator=(const Kernel&) = delete;
  Kernel(Kernel&&) = delete;
  Kernel& operator=(Kernel&&) = delete;

  // Fill the hoisted persistent slots. Must run once before the first run(),
  // and again after anything invalidating the invariants LowerOptions declared.
  //
  // ADDED (reported, Package 4). Lower.h splits the program into a Precompute
  // and a Run stage, and the split is worthless if nothing triggers the first.
  // It is a Kernel entry point rather than something makeKernel does, for the
  // simple reason that it needs the DataTable and makeKernel does not get one.
  // Explicit rather than lazy on first run(): the cost is a full pass over the
  // point set including every hoisted grid lookup, and it belongs in a profile
  // where it is spent -- see the note on precompute() in Interp.h.
  virtual void precompute(const reader::scripting::DataTable& table) = 0;

  // Evaluate every point of the table.
  //
  // CHANGED (reported): this used to promise internal parallelisation. The
  // interpreter does not parallelise -- its tile scratch is per-instance state,
  // so two threads in one run() would share one buffer. Calling run() on one
  // Kernel from several threads is therefore undefined, and a caller wanting
  // parallelism builds one Kernel per thread. An RTC backend may parallelise
  // internally later; that would be a strengthening of this contract, not a
  // change to it.
  //
  // Guarded: running a program with a Precompute stage before precompute() has
  // filled it reads uninitialised persistent slots, which is a silently wrong
  // result rather than a crash, so it is an error rather than a warning.
  virtual void run(const reader::scripting::DataTable& table) = 0;

  // Evaluate a range, with the column BASES supplied per call.
  //
  // ADDED (reported, Package 5). run(table) evaluates a whole table, which is
  // the wrong shape for the consumer that matters: EasiBoundary::query is
  // called once per face and builds a DataTable each time. Measured, that costs
  // 0.4 us per face to build plus 0.4 us to bind, against 0.3 us to evaluate --
  // so today the setup already costs more than the arithmetic, and with a
  // compiled kernel it would be 87% of the total. Compiling is pointless until
  // the setup leaves the per-face loop.
  //
  // So: bind ONCE against a representative table, then call this per face with
  // only the bases moved. Strides, offsets and types stay in the Binding, where
  // they belong -- they are identical for every face.
  //
  // This is what yateto and TensorForge spell PTR_BASED, and it is the right
  // shape here for the same reason: the indirection sits at FACE granularity,
  // so one pointer load amortises over the face's nodal points. An index list
  // would save four bytes per array per face, which at that amortisation is not
  // measurable.
  virtual void run(const KernelArgs& args) = 0;

  [[nodiscard]] virtual BackendKind kind() const = 0;
};

struct BackendOptions {
  BackendKind preferred{BackendKind::Interpreter};
  bool allowFallback{true}; // fall back to Interpreter if `preferred` is unusable
  std::size_t tileSize{0};  // 0 = pick from the program's live-value count
  /// Target architecture for `preferred`, and for that backend ONLY: "native"
  /// or "skylake" for RtcCpu, "80" for NVRTC, "gfx90a" for HIPRTC.
  ///
  /// It is per-backend by nature, so makeKernel CLEARS it when it falls back to
  /// a different one -- otherwise a GPU run that falls back to the compiled CPU
  /// kernel hands "-march=80" to the host compiler, which fails the compile and
  /// silently drops another level to the interpreter. Empty means each backend
  /// picks its own default.
  std::string arch;

  // ADDED (reported, Package 4). makeKernel owns the lowering -- a backend
  // consumes a LoweredProgram, not a Program -- so the only way a caller can
  // declare what is invariant is through here. Lower.h makes hoisting opt-in
  // precisely because the failure mode of getting it wrong is a stale value
  // rather than a crash; leaving these options out of the interface would have
  // made the safe default the only reachable one, and the analytic boundary
  // condition in Package 6 is the case the hoisting exists for.
  LowerOptions lowering;
};

// Chooses a backend, compiles if needed, and returns the prepared kernel.
// Never returns null: on a compilation failure with allowFallback it logs a
// warning and returns the interpreter.
// CHANGED (reported): `binding` is non-const, and `grids` is interned into.
//
// The persistent buffer is sized here because its slot count is a property of
// the LOWERING (declared state plus hoisted values) and Binding::bind() sees
// only the Program. Binding owns the storage -- Program.h puts it next to the
// point set and the permutation it is indexed by -- so sizing it means writing
// to the Binding. Re-sizing is a no-op when the shape already matches, which is
// what lets two kernels over one Binding share the state rather than the second
// one resetting it.
//
// The program's grids are interned into `grids` and the store is loaded here,
// so that a Lookup's dimension and component can be checked against the real
// file once, loudly, instead of being discovered per tile. A caller building
// several kernels over one store should intern all programs first: load() is
// idempotent for opening a file but re-runs the window sizing each time.
std::unique_ptr<Kernel> makeKernel(const Program& program,
                                   Binding& binding,
                                   reader::datafield::GridStore& grids,
                                   const BackendOptions& options);

// --- deferred backends -------------------------------------------------------
// Present as declarations so the decision is recorded in the interface rather
// than in a commit message. Both currently logError on construction.
//
// Texture: hardware trilinear sampling of a device-resident grid. Deferred
// because it is fp32-only and the interpolation weights are quantised to 8–9
// bits of subpixel precision, which breaks CPU/GPU reproducibility — acceptable
// for a material lookup, not acceptable as a default.
//
// Distributed: grids sharded across ranks rather than replicated. Deferred
// because full replication in host memory covers the current models, and
// sharding needs the sub-box/halo analysis on Lookup coordinate arguments
// before it can be correct.
std::unique_ptr<Kernel> makeTextureKernel(const Program& program,
                                          Binding& binding,
                                          reader::datafield::GridStore& grids);
std::unique_ptr<Kernel> makeDistributedKernel(const Program& program,
                                              Binding& binding,
                                              reader::datafield::GridStore& grids);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_BACKEND_H_

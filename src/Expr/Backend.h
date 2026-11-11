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

#include "Expr/Binding.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace seissol::expr {

enum class BackendKind : std::uint8_t {
  Interpreter, // tiled batch/SoA walk; always available, the correctness oracle
  RtcCpu,      // host compiler -> shared object
  RtcCuda,     // NVRTC
  RtcHip,      // hiprtc
  /// OpenCL C source, enqueued through SYCL's kernel_compiler extension.
  ///
  /// The source language is OpenCL C rather than SYCL because that half is
  /// stable; the enqueue goes through SYCL so it works on a Level Zero backend,
  /// where raw OpenCL would not. Support is queried at run time -- whether a
  /// backend accepts OpenCL C source is a property of the backend, not of the
  /// build.
  RtcOpenCl,
  Texture,    // NOT IMPLEMENTED — hardware-sampled grids, see stub below
  Distributed // NOT IMPLEMENTED — sharded grids, see stub below
};

const char* name(BackendKind kind);

/// The backend of the device this build runs on: RtcCpu in a build without one.
[[nodiscard]] constexpr BackendKind deviceBackendOfBuild() {
#if defined(SEISSOL_KERNELS_CUDA)
  return BackendKind::RtcCuda;
#elif defined(SEISSOL_KERNELS_HIP)
  return BackendKind::RtcHip;
#elif defined(SEISSOL_KERNELS_SYCL)
  return BackendKind::RtcOpenCl;
#else
  return BackendKind::RtcCpu;
#endif
}

/// One call's worth of raw arguments, in Binding slot order.
struct KernelArgs {
  /// inputs[i] is the base for Program::inputs()[i], outputs[j] for outputs()[j].
  /// A null entry means "use the base the Binding was bound with", so a caller
  /// that only moves some of the columns need not restate the others.
  const void* const* inputs{nullptr};
  void* const* outputs{nullptr};
  std::size_t inputCount{0};
  std::size_t outputCount{0};

  /// The same for the matrices and blocks of the contractions, in Program::matrices() and
  /// Program::blocks() order. Moving the matrix base is how the subcells of a refined output share
  /// one program: only the pointer changes between the calls.
  const void* const* matrices{nullptr};
  const void* const* blocks{nullptr};
  std::size_t matrixCount{0};
  std::size_t blockCount{0};

  /// The same for the states a table keeps (DataTable::bindState), in Program::state() order. A
  /// state the Binding keeps itself ignores its entry.
  void* const* states{nullptr};
  std::size_t stateCount{0};

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
  // Lower.h splits the program into a Precompute and a Run stage, and the split
  // is worthless if nothing triggers the first. It is a Kernel entry point
  // rather than something makeKernel does, for the simple reason that it needs
  // the DataTable and makeKernel does not get one.
  // Explicit rather than lazy on first run(): the cost is a full pass over the
  // point set including every hoisted grid lookup, and it belongs in a profile
  // where it is spent -- see the note on precompute() in Interp.h.
  virtual void precompute(const reader::scripting::DataTable& table) = 0;

  // Evaluate every point of the table.
  //
  // Not parallelised internally: the tile scratch of the interpreter is
  // per-instance state, so two threads in one run() would share one buffer.
  // Calling run() on one Kernel from several threads is therefore undefined, and
  // a caller wanting parallelism builds one Kernel per thread. An RTC backend may
  // parallelise internally later; that would be a strengthening of this
  // contract, not a change to it.
  //
  // Guarded: running a program with a Precompute stage before precompute() has
  // filled it reads uninitialised persistent slots, which is a silently wrong
  // result rather than a crash, so it is an error rather than a warning.
  virtual void run(const reader::scripting::DataTable& table) = 0;

  // Evaluate a range, with the column BASES supplied per call.
  //
  // run(table) evaluates a whole table, which is the wrong shape for a consumer
  // that asks for a few points many times, such as a boundary condition per face
  // and time step. Building a DataTable per face costs about 0.4 us, binding it
  // another 0.4 us, against 0.3 us to evaluate it -- the setup would cost more
  // than the arithmetic, and with a compiled kernel it would be 87% of the
  // total.
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

/// What the check of a compiled backend found at the first call of its kernel (the specification's
/// check of a backend): the points of the call it compared, their inputs ([input][point], in the
/// order of the binding's inputs), the outputs the interpreter gives and those the backend wrote
/// ([output][point]), how many outputs deviate in their bits, and the largest deviation.
struct BackendComparison {
  std::vector<std::size_t> points;
  std::vector<double> inputs;
  std::vector<double> interpreter;
  std::vector<double> backend;
  std::size_t deviating{0};
  std::string worst;
  /// Compared against the perturbed reference (a backend on a GPU): the tolerance of each value,
  /// [output][point], c = 2 times the largest deviation of R = 8 perturbed evaluations from the
  /// reference; the points where one of them decides a select otherwise than the reference
  /// (unstable); and the largest ratio of a deviation of the backend to its tolerance.
  bool perturbed{false};
  std::vector<double> tolerance;
  std::vector<std::size_t> unstable;
  double worstRatio{0.0};
};

/// Compares what a kernel wrote for `program` in the call of `table` (or of `args`, with `table`
/// null) with what the interpreter gives -- the program as written, nothing hoisted -- at a sample
/// of the points of the call. Nothing for a program with state, which the evaluation would
/// advance, or for outputs that cannot be read back. makeKernel() checks every compiled kernel so
/// at its first call; SEISSOL_EXPR_CHECK=error (the default) stops a run where it deviates, `warn`
/// warns, `off` checks nothing.
///
/// `perturbed` compares as for a backend on a GPU (the specification's "Checking a GPU backend"):
/// against the tolerance of each value from perturbed evaluations, rather than bit for bit.
[[nodiscard]] std::optional<BackendComparison>
    compareWithInterpreter(const Program& program,
                           const Binding& binding,
                           reader::datafield::GridStore& grids,
                           const reader::scripting::DataTable* table,
                           const KernelArgs* args,
                           bool perturbed = false);

/// The same, from the inputs and outputs of the points of a sample read back by whoever can
/// (a backend on a device): inputs[input * points.size() + k] in the order of Program::inputs(),
/// backend[output * points.size() + k]. Nothing for a program with state or contractions.
[[nodiscard]] std::optional<BackendComparison> compareValues(const Program& program,
                                                             reader::datafield::GridStore& grids,
                                                             const std::vector<std::size_t>& points,
                                                             const std::vector<double>& inputs,
                                                             const std::vector<double>& backend,
                                                             bool perturbed);

/// Where the evaluation records of the programs go (the specification's "Evaluation record"): one
/// file per program and backend, written at the first call of its kernel, with the program in the
/// normal form and its fingerprint, the build (`build`, lines as they are), the backend and the
/// compute type, what the check of the backend found, and the sample of points with their inputs
/// and the outputs of the interpreter and of the backend, in hexadecimal. Empty writes none.
void setRecordDirectory(const std::string& directory, const std::vector<std::string>& build = {});

/// Whether the compiled kernels are checked at their first call (SEISSOL_EXPR_CHECK is not off).
[[nodiscard]] bool backendChecks();

/// Logs what the check of the first call of the `kind` kernel of `program` found, and stops the
/// run where it deviates, unless SEISSOL_EXPR_CHECK=warn.
void reportComparison(const Program& program, BackendKind kind, const BackendComparison& found);

/// The points of a call of `count` points from `first` the check compares: all of them up to 64,
/// else 64 spread evenly over the call, its first and its last among them.
[[nodiscard]] std::vector<std::size_t> samplePoints(std::size_t first, std::size_t count);

/// How a compiled CPU kernel takes the points of a tile: one point through the whole program after
/// the other, or one instruction over a block of points after the other (see RtcCpu.h). Auto
/// chooses by the depth of the program.
enum class KernelShape : std::uint8_t { Auto, Point, Tile };

struct BackendOptions {
  BackendKind preferred{BackendKind::Interpreter};
  bool allowFallback{true}; // fall back to Interpreter if `preferred` is unusable
  std::size_t tileSize{0};  // 0 = pick from the program's live-value count
  /// The device queue, for backends that need one to compile as well as to
  /// launch: a `sycl::queue*` for BackendKind::RtcOpenCl. Ignored by every
  /// other backend, which take their stream per call through KernelArgs.
  ///
  /// It is here rather than in KernelArgs because the OpenCL path builds its
  /// kernel bundle from the queue's CONTEXT, which happens once at prepare()
  /// and not per launch.
  void* deviceQueue{nullptr};

  /// Target architecture for `preferred`, and for that backend ONLY: "native"
  /// or "skylake" for RtcCpu, "80" for NVRTC, "gfx90a" for HIPRTC.
  ///
  /// It is per-backend by nature, so makeKernel CLEARS it when it falls back to
  /// a different one -- otherwise a GPU run that falls back to the compiled CPU
  /// kernel hands "-march=80" to the host compiler, which fails the compile and
  /// silently drops another level to the interpreter. Empty means each backend
  /// picks its own default.
  std::string arch;

  // makeKernel owns the lowering -- a backend consumes a LoweredProgram, not a
  // Program -- so the only way a caller can declare what is invariant is through
  // here. Lower.h makes hoisting opt-in precisely because the failure mode of
  // getting it wrong is a stale value rather than a crash; the analytic boundary
  // condition, with its points fixed and its time moving, is the case the
  // hoisting exists for.
  LowerOptions lowering;

  /// Skip the summary line a kernel logs when it is made: for the second and later kernels of
  /// one program, e.g. one per thread.
  bool quiet{false};

  /// The shape of a compiled CPU kernel; for the RtcCpu backend only.
  KernelShape shape{KernelShape::Auto};
};

// Chooses a backend, compiles if needed, and returns the prepared kernel.
// Never returns null: on a compilation failure with allowFallback it logs a
// warning and returns the interpreter. Writes to `binding` and `grids`, see
// below.
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
// file once, loudly, instead of being discovered per tile. A caller that makes
// its kernels on some ranks only interns the grids of its program on all ranks
// before (internGrids()): loading is collective where it opens or sizes a grid
// (GridStore::load()).
std::unique_ptr<Kernel> makeKernel(const Program& program,
                                   Binding& binding,
                                   reader::datafield::GridStore& grids,
                                   const BackendOptions& options);

/// Interns the grids of `program` into `grids` and loads them; their indices in the store, by
/// GridId. Collective where it opens or sizes a grid (GridStore::load()): a script whose kernels
/// are made on some ranks only -- where a cluster has a face of its boundary, say -- calls this
/// on all ranks when it is loaded, so that the store opens the same grids in the same order on
/// all of them, and making a kernel later finds them open.
std::vector<std::size_t> internGrids(const Program& program, reader::datafield::GridStore& grids);

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

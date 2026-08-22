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

// A prepared, ready-to-run instance of one Program against one Binding.
class Kernel {
  public:
  Kernel() = default;
  virtual ~Kernel() = default;
  Kernel(const Kernel&) = delete;
  Kernel& operator=(const Kernel&) = delete;
  Kernel(Kernel&&) = delete;
  Kernel& operator=(Kernel&&) = delete;

  // Evaluate every point of the table. Thread-safety: `run` may be called from
  // one thread only; it parallelises internally.
  virtual void run(const reader::scripting::DataTable& table) = 0;

  [[nodiscard]] virtual BackendKind kind() const = 0;
};

struct BackendOptions {
  BackendKind preferred{BackendKind::Interpreter};
  bool allowFallback{true}; // fall back to Interpreter if `preferred` is unusable
  std::size_t tileSize{0};  // 0 = pick from the program's live-value count
  std::string arch;         // "native", "sm_80", "gfx90a", …
};

// Chooses a backend, compiles if needed, and returns the prepared kernel.
// Never returns null: on a compilation failure with allowFallback it logs a
// warning and returns the interpreter.
std::unique_ptr<Kernel> makeKernel(const Program& program,
                                   const Binding& binding,
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
                                          const Binding& binding,
                                          reader::datafield::GridStore& grids);
std::unique_ptr<Kernel> makeDistributedKernel(const Program& program,
                                              const Binding& binding,
                                              reader::datafield::GridStore& grids);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_BACKEND_H_

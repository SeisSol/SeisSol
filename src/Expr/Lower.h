// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_LOWER_H_
#define SEISSOL_SRC_EXPR_LOWER_H_

// DAG -> straight-line SSA, with liveness-based slot reuse. This is the pass the
// handover calls "eigenes Lowering"; it is NOT the sderiv lower.hpp, which is
// the derived-output cadence/region planner and stays where it is.
//
// Two things make this more than a flattening:
//
// SLOT REUSE. A tile is meant to sit in L1. Without reuse, ~30 nodes at f64 and
// 256 lanes is 60 KiB and the tile is a DRAM structure with extra steps. Values
// die at their last use, so the live count is what actually has to fit, and the
// backend picks the tile width from it rather than the other way round.
//
// STAGES. Values that cannot change between run() calls are computed once, in
// prepare(), and read back from persistent storage afterwards — loop-invariant
// code motion across the timestep loop rather than inside it. For the analytic
// boundary condition this is the whole ballgame: x, y and z are fixed for the
// life of the run, so every grid Lookup on them is invariant, and the width^d
// random-access gather that Paket 3 identifies as the expensive part happens
// once instead of once per timestep.
//
// Hoisting is opt-in and conservative by default: LowerOptions names what the
// caller *knows* to be constant, and an empty option set hoists nothing. The
// failure mode of getting this wrong is a stale value, not a crash, so the
// default has to be the safe one.

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace seissol::expr {

// Index space a value lives in. Pointwise ops require their operands in one
// space; Uniform is the bottom element and broadcasts into any of them.
//
// Today there is exactly one concrete space, so the join below can never fail.
// It is computed and carried anyway, because it is the hook a nodal<->modal
// transform needs: such a node is a reduction over the batch dimension, i.e.
// every output lane reads every input lane, which makes it a tile boundary by
// construction and gives the regions on either side of it their own extent.
// Adding it later is additive if the tag is already threaded through; it is a
// rewrite of the allocator if it is not.
enum class Space : std::uint8_t { Uniform, Point };

const char* name(Space space);

// Lifetime class. Transient values live inside one tile and share slots;
// persistent values are point-indexed and survive across run() calls. The two
// have different owners: transients belong to the backend's tile buffer,
// persistents to Binding, which is the layer that knows numPoints() and owns the
// permutation the indices are relative to.
enum class ValueClass : std::uint8_t { Transient, Persistent };

enum class Opcode : std::uint8_t {
  Const,          // materialise Instruction::value
  LoadInput,      // imm = input index in the gathered tile
  LoadPersistent, // imm = persistent slot (state slot or hoisted value)
  Pw,             // fn, arity(fn) operands
  Lookup          // grid, comp, argCount coordinate operands
};

const char* name(Opcode op);

// Every operand refers to a TRANSIENT slot of the stage it appears in. Values
// that cross a stage boundary do so through an explicit LoadPersistent /
// persistent Store pair rather than by letting operands address two different
// memories. That costs one copy per crossing value per tile — negligible next to
// the memory traffic the load itself causes — and buys a uniform inner loop and
// a liveness problem that stays straight-line.
struct Instruction {
  Opcode op{Opcode::Const};
  Fn fn{Fn::Neg};
  Space space{Space::Point};
  std::int32_t imm{0};
  double value{0.0};
  GridId grid{NoGrid};
  std::int32_t comp{0};
  std::int32_t operandBegin{0};
  std::int32_t operandCount{0};
  std::int32_t dst{0}; // transient slot
  NodeId node{NoNode}; // provenance, for diagnostics and codegen comments
};

struct Store {
  std::int32_t target{0}; // output index, or persistent slot
  std::int32_t source{0}; // transient slot
  Space space{Space::Point};
};

struct StageCode {
  std::vector<Instruction> code;
  std::vector<Store> outputs;    // Run only: scatter targets, by output index
  std::vector<Store> persistent; // Precompute: hoisted values. Run: state writeback.
  std::int32_t slotCount{0};     // transient slots this stage needs
};

inline constexpr std::int32_t DefaultHoistThreshold = 4;

struct LowerOptions {
  // Inputs whose value is identical on every run() call for one Binding.
  std::vector<std::string> invariantInputs;
  // Grids not touched by GridStore::update(). A grid left out of this list is
  // never hoisted, which is the safe direction: hoisting an updating grid is
  // silently stale, not loud.
  std::vector<GridId> invariantGrids;
  // Minimum estimated cost before an invariant value is worth materialising
  // instead of recomputing. Materialising costs numPoints * sizeof(ComputeType)
  // of persistent memory and a streaming read per call, so recomputation wins
  // for cheap nodes. The weights behind this are a placeholder for a measurement.
  std::int32_t hoistThreshold{DefaultHoistThreshold};

  // Identity of the lowering configuration. NOT part of Program::fingerprint():
  // options shape the kernel but not the program, so the kernel cache mixes this
  // in next to arch and backend, the way cache.hpp already mixes those.
  [[nodiscard]] std::uint64_t fingerprint() const;
};

class LoweredProgram {
  public:
  [[nodiscard]] const StageCode& precompute() const { return precompute_; }
  [[nodiscard]] const StageCode& run() const { return run_; }
  [[nodiscard]] const std::vector<std::int32_t>& operands() const { return operands_; }

  // Persistent slots [0, stateSlotCount) are Program::state(), in declaration
  // order, and are the ones Binding initialises from StateSpec::initial. The
  // remainder are hoisted values, which precompute() writes before any run().
  [[nodiscard]] std::int32_t stateSlotCount() const { return stateSlotCount_; }
  [[nodiscard]] std::int32_t persistentSlotCount() const { return persistentSlotCount_; }
  [[nodiscard]] bool hasPrecompute() const { return !precompute_.code.empty(); }

  // Peak transient slots over both stages — the number the backend divides its
  // L1 budget by to get the tile width.
  [[nodiscard]] std::int32_t peakSlotCount() const;

  // One line for prepare() to log. The difference between a hoisted and a
  // non-hoisted kernel is large enough that it must not be invisible.
  [[nodiscard]] std::string summary() const;
  [[nodiscard]] std::string dump() const;

  private:
  friend LoweredProgram lower(const Program&, const LowerOptions&);
  StageCode precompute_;
  StageCode run_;
  std::vector<std::int32_t> operands_;
  std::int32_t stateSlotCount_{0};
  std::int32_t persistentSlotCount_{0};
};

// Throws std::invalid_argument on a program that validate() would reject and on
// an operand space mismatch. Call validate() first; this does not repeat it.
LoweredProgram lower(const Program& program, const LowerOptions& options = {});

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_LOWER_H_

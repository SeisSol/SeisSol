// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_LOWER_H_
#define SEISSOL_SRC_EXPR_LOWER_H_

// DAG -> straight-line SSA, with liveness-based slot reuse.
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
// random-access gather, which is the expensive part, happens once instead of
// once per timestep.
//
// Hoisting is opt-in and conservative by default: LowerOptions names what the
// caller *knows* to be constant, and an empty option set hoists nothing. The
// failure mode of getting this wrong is a stale value, not a crash, so the
// default has to be the safe one.
//
// BRANCHES. The values that only the points taking one branch of a select need
// come together in a part of the code a backend may skip where no point of its
// tile, warp or own takes the branch (StageCode::branches, Branches.h). The code
// stays straight-line and valid as it is: a backend that runs every part
// computes the same values, and a select reads the operand it does not choose
// from a slot that holds some value of its type.

#include "Expr/Expand.h"
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
  GridParam,      // grid, comp = GridField, imm = axis; no operands (Kind::GridParam)
  GridLoad,       // grid; one operand, the offset (Kind::GridLoad)
  Contract,       // matrix, block; no operands
  Convert,        // one operand, converted from the type of its slot to the type of the instruction
  Percentile      // the percentile, then the values (Kind::Percentile)
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
  // The type the instruction computes in, and its destination holds. A load converts into it from
  // the type the backend moves the data in (LoweredProgram::dataType()); the operands of a Pw, a
  // Lookup and a Contract have it as well (typed() sees to it), the operand of a Convert has
  // another.
  ComputeType type{ComputeType::F64};
  Fn fn{Fn::Neg};
  Space space{Space::Point};
  std::int32_t imm{0};
  double value{0.0};
  GridId grid{NoGrid};
  std::int32_t comp{0};
  MatrixId matrix{NoMatrix};
  BlockId block{NoBlock};
  std::int32_t operandBegin{0};
  std::int32_t operandCount{0};
  std::int32_t dst{0}; // transient slot
  NodeId node{NoNode}; // provenance, for diagnostics and codegen comments
  // Pw: computed as written, without contraction (Node::strict)
  bool strict{false};
};

// A write of a transient into an output or a persistent slot. It runs before the instruction at
// `position` (at the end of the code where that is its size), outside every branch, so at every
// point: as soon as its value is computed, and the value it would overwrite loaded -- the old value
// of a state, an input of the same name as an output --, so that a value that is written lives no
// longer than needed (see THE ORDER in Branches.h). The writes of a position run in the order of
// their lists, outputs first.
struct Store {
  std::int32_t target{0}; // output index, or persistent slot
  std::int32_t source{0}; // transient slot
  Space space{Space::Point};
  std::int32_t position{0};
};

/// A part of the code of a stage that only the points where the condition is `side` need: the
/// instructions [begin, end), the values of the branch of a select (Branches.h). A backend runs it
/// for a set of points where one of them needs it, and may skip it where none does; it reads the
/// condition before the part. The parts nest: one that starts where another does, or in it, ends
/// in it.
struct Branch {
  std::int32_t condition{0}; // transient slot, of any type: true where not zero, a NaN included
  bool side{true};
  std::int32_t begin{0};
  std::int32_t end{0};
};

struct CheckStore {
  GridId grid{NoGrid};
  CheckKind kind{CheckKind::Outside};
  std::int32_t source{0}; // transient slot
  std::int32_t index{-1}; // of kind User: the check of the program
};

struct StageCode {
  std::vector<Instruction> code;
  // Run only: scatter targets, by output index
  std::vector<Store> outputs;
  // Precompute: hoisted values. Run: state writeback, none for a state whose next value is its
  // value.
  std::vector<Store> persistent;
  // Run only: the checks of the written-out lookups (GridCheck), each of a grid and a kind, read
  // from a transient of type Bool that is true where the check finds what its kind says; after the
  // code.
  std::vector<CheckStore> checks;
  // By `begin`; of two that start at one instruction, the one that holds the other first.
  std::vector<Branch> branches;
  std::int32_t slotCount{0}; // transient slots this stage needs
  // The type of each transient slot: a slot only ever holds values of one type. A store converts
  // from it into the type the data is moved in.
  std::vector<ComputeType> slotTypes;
};

inline constexpr std::int32_t DefaultHoistThreshold = 4;
// About four elements of a grid, or two operations of the maths library: what a branch skipped
// for a whole tile or warp saves, against the test of its condition, and a GPU warp whose points
// diverge, which runs both.
inline constexpr std::int32_t DefaultBranchThreshold = 16;

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
  // Minimum weight, by the same estimate, of the values that only one branch of a select needs,
  // before they become a branch of the code (StageCode::branches). A larger value than any program
  // weighs makes no branches, which gives the same values.
  std::int32_t branchThreshold{DefaultBranchThreshold};

  // Identity of the lowering configuration. NOT part of Program::fingerprint():
  // options shape the kernel but not the program, so the kernel cache mixes this
  // in next to arch and backend.
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

  /// The compute type of the program, and the type the backends move its data in: the inputs and
  /// outputs, the states and the hoisted values (see dataType()).
  [[nodiscard]] ComputeType computeType() const { return computeType_; }
  [[nodiscard]] ComputeType dataType() const { return dataType_; }

  /// Whether an instruction computes in another type than the compute type.
  [[nodiscard]] bool mixesTypes() const;
  [[nodiscard]] std::int32_t persistentSlotCount() const { return persistentSlotCount_; }
  [[nodiscard]] bool hasPrecompute() const { return !precompute_.code.empty(); }

  // Peak transient slots over both stages — the number the backend divides its
  // L1 budget by to get the tile width.
  [[nodiscard]] std::int32_t peakSlotCount() const;

  // One line for prepare() to log. The difference between a hoisted and a
  // non-hoisted kernel is large enough that it must not be invisible.
  [[nodiscard]] std::string summary() const;
  [[nodiscard]] std::string dump() const;

  /// How the program samples each of its grids (by GridId), and whether it samples any: the
  /// lookups are written out (expandLookups()), so this is what a backend checks the grids it is
  /// bound to against.
  [[nodiscard]] const std::vector<GridUse>& gridUses() const { return gridUses_; }
  [[nodiscard]] bool usesGrids() const;
  /// Whether an instruction is computed without contraction (Instruction::strict).
  [[nodiscard]] bool hasStrict() const;
  /// Whether an instruction is a percentile, which a compiled backend computes with a function of
  /// its own (Opcode::Percentile).
  [[nodiscard]] bool usesPercentile() const;

  private:
  friend LoweredProgram lower(const Program&, const LowerOptions&);
  StageCode precompute_;
  StageCode run_;
  std::vector<std::int32_t> operands_;
  std::int32_t stateSlotCount_{0};
  std::int32_t persistentSlotCount_{0};
  ComputeType computeType_{ComputeType::F64};
  ComputeType dataType_{ComputeType::F64};
  std::vector<GridUse> gridUses_;
};

// Throws std::invalid_argument on a program that validate() would reject and on
// an operand space mismatch. Call validate() first; this does not repeat it. The program is lowered
// with the type of every value resolved (typed()) and its lookups written out (expandLookups());
// its inputs, outputs and states stay in order.
LoweredProgram lower(const Program& program, const LowerOptions& options = {});

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_LOWER_H_

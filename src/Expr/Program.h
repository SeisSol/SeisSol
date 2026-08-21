// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_PROGRAM_H_
#define SEISSOL_SRC_EXPR_PROGRAM_H_

// A `Program` is what a frontend produces and a backend consumes: an arena, one
// DAG root per named output, the inferred input/output signature, and the list
// of external data grids the DAG references.
//
// SIGNATURE IS INFERRED, NOT DECLARED. The current readers take
// `M.input_parameters` / `M.output_parameters` on trust and match them against
// DataTable column names inside every call(); every mismatch found so far in
// that code was a variant of "the name silently resolved to the wrong column".
// Here the frontend reports what the DAG actually reads and writes, and Binding
// checks it against the table once.
//
// TYPING: there is deliberately no per-node type. Everything inside a Program is
// computed in ONE type, chosen per Program; conversion happens only at the
// DataTable boundary (getValueAs/setValueAs). The reason not to build a real
// type system is that the only mixed-type inputs in practice are `real` columns
// (f32 in single-precision builds) and the integer `group`/`sim` columns, and
// both are handled correctly by a single compute type plus exact conversion at
// the edges. The one caveat is worth writing down rather than solving: with
// ComputeType::F32, integer channels above 2^24 lose exactness. Group and fault
// tag ids are far below that; if that ever stops being true, this is the
// assumption that breaks, not something subtler.

#include "Expr/Ir.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace seissol::expr {

// The type every intermediate value is computed in. F32 exists for the
// per-timestep GPU path, where fp64 throughput is 1/32–1/64 on most parts.
enum class ComputeType : std::uint8_t { F32, F64 };

struct VarSpec {
  std::string name;
  reader::scripting::DataType type{reader::scripting::DataType::F64};
};

// ADDED (reported). A state slot is a program-owned, point-indexed value that
// survives across run() calls. It exists because a script can *invent* an
// accumulator that no consumer knows about: `acc = acc + f(x,t)*dt` has no
// DataTable column to bind to, and Binding::bind would reject it as "required
// input with no matching column". That would break exactly the property this
// header is built around — the signature being inferred rather than declared.
//
// This is NOT Kind::Cumint by another name. The IR stays pure: the state is read
// as an ordinary Field channel and written by an ordinary DAG root. What is
// stateful is the storage, which lives in Binding, next to the point set whose
// identity it depends on, and next to the permutation that reorders it.
//
// SEMANTICS ARE PARALLEL ASSIGNMENT, NOT SEQUENTIAL. Every state root is
// evaluated against the values from the *previous* call, then all state slots
// are written. Two states referring to each other are therefore well defined and
// order-independent. This falls out of the lowering (all stores happen at the
// end of a tile), but it is a promise, not an accident.
struct StateSpec {
  std::string name;    // the Field channel name the DAG reads it back through
  double initial{0.0}; // value on the first call, and after every rebind
  NodeId root{NoNode}; // DAG root producing the value for the next call
};

class Program {
  public:
  [[nodiscard]] const Arena& arena() const { return arena_; }
  [[nodiscard]] Arena& arena() { return arena_; }

  [[nodiscard]] const std::vector<VarSpec>& inputs() const { return inputs_; }
  [[nodiscard]] const std::vector<VarSpec>& outputs() const { return outputs_; }
  [[nodiscard]] const std::vector<StateSpec>& state() const { return state_; }

  // roots()[i] is the DAG root producing outputs()[i].
  [[nodiscard]] const std::vector<NodeId>& roots() const { return roots_; }

  // Every external grid the DAG references, in GridId order.
  [[nodiscard]] const std::vector<datafield::GridDesc>& grids() const { return grids_; }

  [[nodiscard]] ComputeType computeType() const { return computeType_; }
  void setComputeType(ComputeType t) { computeType_ = t; }

  // Canonical byte form of everything the emitted kernel depends on, and
  // nothing else. Deterministic across runs and builds: the walk is a
  // depth-first post-order from the roots in root order, visiting children in
  // Arena::children() order, so the result is a function of the DAG's *shape*
  // and not of NodeId values, which in turn depend on frontend construction
  // order. Channel names are absent — Binding resolves them to slot indices, so
  // they do not reach the emitter; a Field is serialised as its index in
  // inputs() or state(). Grid file paths are absent for the same reason, which
  // is why GridDesc has to supply canonicalKey().
  [[nodiscard]] std::string canonicalForm() const;

  // Stable, backend-independent fingerprint of the whole program. Used as the
  // kernel-cache key and as the identity for "have I already compiled this".
  //
  // CHANGED (reported): this is the hash of canonicalForm(), and the cache is
  // expected to keep the bytes and compare them on a hit rather than trusting
  // 64 bits. A fingerprint collision does not produce a diagnostic, it produces
  // a silently wrong kernel; the string compare costs nothing at cache-lookup
  // frequency. Lowering options (LICM configuration) are deliberately NOT in
  // here — they shape the kernel but not the program, and the cache mixes them
  // in alongside arch and backend, exactly as cache.hpp already does.
  [[nodiscard]] std::uint64_t fingerprint() const;

  // --- construction (frontends only) ---
  void addOutput(const std::string& name, reader::scripting::DataType type, NodeId root);
  void addInput(const std::string& name, reader::scripting::DataType type);
  void addState(const std::string& name, double initial, NodeId root);
  GridId internGrid(const datafield::GridDesc& desc);

  private:
  Arena arena_;
  std::vector<VarSpec> inputs_;
  std::vector<VarSpec> outputs_;
  std::vector<StateSpec> state_;
  std::vector<NodeId> roots_;
  std::vector<datafield::GridDesc> grids_;
  ComputeType computeType_{ComputeType::F64};
};

// Structural checks that do not need a DataTable: every Field channel is
// declared as an input or a state, every root is in range, no Lookup references
// an out-of-range GridId, arities match, and the program does not require state
// or an element context that the scripting path cannot supply.
//
// CHANGED (reported): the last clause used to be spelled as a list of forbidden
// kinds (Dx, Cumint, Fold, Sample). It is now expressed through two predicates
// over Kind — requiresImplicitState() and requiresElementContext() — so that a
// new node kind is classified once, at the enum, instead of every site that
// happens to enumerate the forbidden set. Cumint and Fold are rejected because
// they need per-point state and a dt that does not exist at a boundary node
// under LTS; Dx and Sample because they need modal DOFs and an output cadence
// respectively. Declared state (Program::state()) is the supported way to keep
// a value across calls.
//
// Throws std::invalid_argument; the caller turns that into logError.
void validate(const Program& program);

// True for kinds that carry state the scripting ABI has nowhere to put, and for
// kinds that only mean something inside an element/output-cadence context.
[[nodiscard]] constexpr bool requiresImplicitState(Kind kind) {
  return kind == Kind::Cumint || kind == Kind::Fold;
}
[[nodiscard]] constexpr bool requiresElementContext(Kind kind) {
  return kind == Kind::Dx || kind == Kind::Sample;
}

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_PROGRAM_H_

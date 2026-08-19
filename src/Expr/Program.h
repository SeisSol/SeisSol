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

class Program {
  public:
  [[nodiscard]] const Arena& arena() const { return arena_; }
  [[nodiscard]] Arena& arena() { return arena_; }

  [[nodiscard]] const std::vector<VarSpec>& inputs() const { return inputs_; }
  [[nodiscard]] const std::vector<VarSpec>& outputs() const { return outputs_; }

  // roots()[i] is the DAG root producing outputs()[i].
  [[nodiscard]] const std::vector<NodeId>& roots() const { return roots_; }

  // Every external grid the DAG references, in GridId order.
  [[nodiscard]] const std::vector<datafield::GridDesc>& grids() const { return grids_; }

  [[nodiscard]] ComputeType computeType() const { return computeType_; }
  void setComputeType(ComputeType t) { computeType_ = t; }

  // Stable, backend-independent fingerprint of the whole program. Used as the
  // kernel-cache key and as the identity for "have I already compiled this".
  [[nodiscard]] std::uint64_t fingerprint() const;

  // --- construction (frontends only) ---
  void addOutput(const std::string& name, reader::scripting::DataType type, NodeId root);
  void addInput(const std::string& name, reader::scripting::DataType type);
  GridId internGrid(const datafield::GridDesc& desc);

  private:
  Arena arena_;
  std::vector<VarSpec> inputs_;
  std::vector<VarSpec> outputs_;
  std::vector<NodeId> roots_;
  std::vector<datafield::GridDesc> grids_;
  ComputeType computeType_{ComputeType::F64};
};

// Structural checks that do not need a DataTable: every Field channel is
// declared as an input, every root is reachable, no Lookup references an
// out-of-range GridId, arities match, no derived-output-only node (Dx, Cumint,
// Fold, Sample) appears in a Program destined for the scripting path.
// Throws std::invalid_argument; the caller turns that into logError.
void validate(const Program& program);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_PROGRAM_H_

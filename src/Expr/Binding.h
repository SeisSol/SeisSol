// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_BINDING_H_
#define SEISSOL_SRC_EXPR_BINDING_H_

// Binding resolves a Program's signature against a concrete DataTable ONCE, in
// DataReader::prepare(), and produces the gather/scatter descriptors the
// backends work from.
//
// This is the layer that removes the per-point cost of the current readers. In
// EasiReader/LuaReader the name lookup, the direction check and the datatype
// dispatch all happen inside the point loop, behind a std::function. Here the
// lookup happens once, the datatype dispatch collapses into a converter chosen
// once per column, and the backend only ever sees contiguous SoA tiles.
//
// GROUP PARTITIONING lives here too, and it is the reason the IR needs no
// dedicated switch/groups node. Select alone would be wrong for layered models:
// a select chain evaluates every branch at every point, so a model with one
// ASAGI grid per layer would read every grid at every point. Instead the tile
// builder partitions the point set by the `group` column and runs each
// sub-batch through the branch that its group selects — the same thing easi's
// Switch does with subsetAdapter, one level lower. Programs without a group
// input get a single partition covering all points.

#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <vector>

namespace seissol::expr {

// How one DataTable column maps into (or out of) a compute-type tile buffer.
struct ColumnBinding {
  std::size_t entry{0}; // index into DataTable::dataEntries()
  int slot{0};          // input channel id / output index
  reader::scripting::DataType tableType{reader::scripting::DataType::F64};
  bool computed{false}; // bindComputed column: gather goes through the functor
};

// A contiguous run of points sharing one group value. Half-open [begin, end).
struct GroupRange {
  std::int32_t group{0};
  std::size_t begin{0};
  std::size_t end{0};
};

class Binding {
  public:
  // Validates the Program against the table and resolves every column.
  // Throws std::invalid_argument on: a required input with no matching column,
  // an output bound to an In-only column, a duplicate column name, or a point
  // count of zero. Note that an *extra* table column is not an error — the
  // consumer is allowed to offer more than the program reads.
  static Binding bind(const Program& program, const reader::scripting::DataTable& table);

  [[nodiscard]] const std::vector<ColumnBinding>& inputs() const { return inputs_; }
  [[nodiscard]] const std::vector<ColumnBinding>& outputs() const { return outputs_; }
  [[nodiscard]] std::size_t numPoints() const { return numPoints_; }

  // Present only when the program reads a `group` channel; empty otherwise.
  [[nodiscard]] const std::vector<GroupRange>& groupRanges() const { return groupRanges_; }
  // Permutation applied to reach the group ranges; empty when unpermuted.
  [[nodiscard]] const std::vector<std::size_t>& permutation() const { return permutation_; }

  // Gather `count` points starting at `first` into `dst`, one contiguous lane
  // block per input channel: dst[channel * count + lane]. Scatter is the
  // mirror image. Both are the only places a DataTable accessor is touched.
  //
  // Integer channels are NOT overloaded, and that is deliberate. A group or
  // fault-tag column is read through the same tile as everything else, because
  // the program has one compute type and Program.h already records what that
  // costs (exactness above 2^24 under F32). Adding int32/int64 tiles would
  // reintroduce the per-node typing that Program.h declines to build.
  void gather(const reader::scripting::DataTable& table,
              std::size_t first,
              std::size_t count,
              double* dst) const;
  void gather(const reader::scripting::DataTable& table,
              std::size_t first,
              std::size_t count,
              float* dst) const;
  void scatter(const reader::scripting::DataTable& table,
               std::size_t first,
               std::size_t count,
               const double* src) const;
  void scatter(const reader::scripting::DataTable& table,
               std::size_t first,
               std::size_t count,
               const float* src) const;

  private:
  std::vector<ColumnBinding> inputs_;
  std::vector<ColumnBinding> outputs_;
  std::vector<GroupRange> groupRanges_;
  std::vector<std::size_t> permutation_;
  std::size_t numPoints_{0};
};

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_BINDING_H_

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
#include <cstdint>
#include <optional>
#include <vector>

namespace seissol::expr {

// The input channel a program reads its layer/region id through. Spelled once,
// here, because Binding is the only layer that treats a channel as anything
// other than a name: it is the column the tile builder partitions on.
inline constexpr const char* GroupChannelName = "group";

// How one DataTable column maps into (or out of) a compute-type tile buffer.
struct ColumnBinding {
  std::size_t entry{0}; // index into DataTable::dataEntries()
  int slot{0};          // index into Program::inputs() / Program::outputs()
  reader::scripting::DataType tableType{reader::scripting::DataType::F64};
  // The column has no setter, so it can never serve as an output. True for
  // bindComputed and for the bindViewConst/bindMemberViewConst family alike --
  // DataTable does not distinguish them from the outside, and the property that
  // matters downstream is the writability, not which builder produced it.
  bool computed{false};

  /// Address arithmetic, resolved once at bind time. Empty exactly when the
  /// column was bound with bindComputed -- see DataTable::StridedView. A
  /// program with any empty entry can only be evaluated through the accessor,
  /// which rules out the device backends and the per-call base override.
  std::optional<reader::scripting::StridedView> view;
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

  // --- persistent storage ---
  //
  // ADDED (reported). Program.h puts the state next to the point set and the
  // permutation it is indexed by, i.e. here; Interp.h expects a raw pointer with
  // slot-major layout, persistent[slot * numPoints + point]. The buffer cannot
  // be sized in bind(), because the slot count is a property of the LOWERING
  // (state slots plus hoisted values) and bind() sees only the Program. Hence a
  // second, explicit call once the lowering exists, rather than a hidden resize
  // on first use: the allocation is numPoints * slots * sizeof(ComputeType) and
  // belongs where a profile can see it.
  //
  // Re-allocating resets the state to StateSpec::initial. That is the documented
  // meaning of a rebind -- a state slot is tied to the identity of the point set
  // it was allocated for, and a new point set has no history to carry.
  // A no-op when the shape already matches, so two kernels over one Binding
  // share the state instead of the second resetting it to StateSpec::initial.
  void allocatePersistent(const Program& program, std::int32_t slotCount);
  [[nodiscard]] std::int32_t persistentSlotCount() const { return persistentSlotCount_; }

  // Typed views of the buffer above. Not a template, for the reason the gather
  // overloads are not: the compute type is fixed per Program, so the choice is
  // made once by the caller that already switched on it to pick an interpreter.
  // Both log an error when asked for the type the Program does not compute in.
  [[nodiscard]] double* persistentF64();
  [[nodiscard]] float* persistentF32();

  // Gather `count` points starting at `first` into `dst`, one contiguous lane
  // block per input channel: dst[channel * count + lane]. Scatter is the
  // mirror image. Both are the only places a DataTable accessor is touched.
  //
  // Integer channels are NOT overloaded, and that is deliberate. A group or
  // fault-tag column is read through the same tile as everything else, because
  // the program has one compute type and Program.h already records what that
  // costs (exactness above 2^24 under F32). Adding int32/int64 tiles would
  // reintroduce the per-node typing that Program.h declines to build.
  /// True when every bound column has a StridedView, i.e. when this binding can
  /// be evaluated from raw pointers alone. What makeKernel checks before
  /// offering a device backend, and what run(KernelArgs) requires.
  [[nodiscard]] bool addressable() const { return addressable_; }

  /// Gather from bases supplied per call rather than from the bound table.
  /// `inputs[i]` may be null to keep the bound base for that slot.
  void gatherFrom(const void* const* inputs,
                  std::size_t inputCount,
                  std::size_t first,
                  std::size_t count,
                  double* dst) const;
  void gatherFrom(const void* const* inputs,
                  std::size_t inputCount,
                  std::size_t first,
                  std::size_t count,
                  float* dst) const;
  void scatterTo(void* const* outputs,
                 std::size_t outputCount,
                 std::size_t first,
                 std::size_t count,
                 const double* src) const;
  void scatterTo(void* const* outputs,
                 std::size_t outputCount,
                 std::size_t first,
                 std::size_t count,
                 const float* src) const;

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
  void buildGroupRanges(const Program& program, const reader::scripting::DataTable& table);

  template <typename Tile>
  void gatherFromImpl(const void* const* inputs,
                      std::size_t inputCount,
                      std::size_t first,
                      std::size_t count,
                      Tile* dst) const;
  template <typename Tile>
  void scatterToImpl(void* const* outputs,
                     std::size_t outputCount,
                     std::size_t first,
                     std::size_t count,
                     const Tile* src) const;

  std::vector<ColumnBinding> inputs_;
  std::vector<ColumnBinding> outputs_;
  std::vector<GroupRange> groupRanges_;
  std::vector<std::size_t> permutation_;
  std::size_t numPoints_{0};
  bool addressable_{false};
  std::vector<std::byte> persistent_;
  std::int32_t persistentSlotCount_{0};
  ComputeType computeType_{ComputeType::F64};
};

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_BINDING_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Expr/Binding.h"

#include "Expr/Interp.h"
#include "Expr/Program.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace seissol::expr {

namespace {

using reader::scripting::DataEntry;
using reader::scripting::DataTable;
using reader::scripting::DataType;
using reader::scripting::Direction;

const char* name(DataType type) {
  switch (type) {
  case DataType::F32:
    return "f32";
  case DataType::F64:
    return "f64";
  case DataType::I32:
    return "i32";
  case DataType::I64:
    return "i64";
  }
  return "unknown";
}

// One column, one lane block. The switch on the column's storage type is OUT of
// the lane loop on purpose: that dispatch is the per-point cost Binding.h exists
// to remove, and DataEntry::getValueAs would put it back by taking it inside.
//
// getValue<Col> is used rather than getValueAs<Tile> for a second reason:
// getValueAs ends in a bare `throw;` outside any catch, which terminates rather
// than diagnoses if the enum ever grows a value. Naming the column type here
// means the exhaustive switch lives at one site instead of at every call.
template <typename Tile, typename Col>
void gatherColumn(const DataEntry& entry,
                  const std::vector<std::size_t>& permutation,
                  std::size_t first,
                  std::size_t count,
                  Tile* dst) {
  if (permutation.empty()) {
    for (std::size_t lane = 0; lane < count; ++lane) {
      dst[lane] = static_cast<Tile>(entry.getValue<Col>(first + lane));
    }
  } else {
    for (std::size_t lane = 0; lane < count; ++lane) {
      dst[lane] = static_cast<Tile>(entry.getValue<Col>(permutation[first + lane]));
    }
  }
}

template <typename Tile>
void gatherOne(const ColumnBinding& column,
               const DataTable& table,
               const std::vector<std::size_t>& permutation,
               std::size_t first,
               std::size_t count,
               Tile* dst) {
  const DataEntry& entry = table.dataEntries()[column.entry];
  switch (column.tableType) {
  case DataType::F32:
    gatherColumn<Tile, float>(entry, permutation, first, count, dst);
    return;
  case DataType::F64:
    gatherColumn<Tile, double>(entry, permutation, first, count, dst);
    return;
  case DataType::I32:
    gatherColumn<Tile, std::int32_t>(entry, permutation, first, count, dst);
    return;
  case DataType::I64:
    gatherColumn<Tile, std::int64_t>(entry, permutation, first, count, dst);
    return;
  }
  logError() << "expr: column" << entry.name << "has an unhandled storage type.";
}

template <typename Tile, typename Col>
void scatterColumn(const DataEntry& entry,
                   const std::vector<std::size_t>& permutation,
                   std::size_t first,
                   std::size_t count,
                   const Tile* src) {
  if (permutation.empty()) {
    for (std::size_t lane = 0; lane < count; ++lane) {
      entry.setValue<Col>(first + lane, static_cast<Col>(src[lane]));
    }
  } else {
    for (std::size_t lane = 0; lane < count; ++lane) {
      entry.setValue<Col>(permutation[first + lane], static_cast<Col>(src[lane]));
    }
  }
}

template <typename Tile>
void scatterOne(const ColumnBinding& column,
                const DataTable& table,
                const std::vector<std::size_t>& permutation,
                std::size_t first,
                std::size_t count,
                const Tile* src) {
  const DataEntry& entry = table.dataEntries()[column.entry];
  switch (column.tableType) {
  case DataType::F32:
    scatterColumn<Tile, float>(entry, permutation, first, count, src);
    return;
  case DataType::F64:
    scatterColumn<Tile, double>(entry, permutation, first, count, src);
    return;
  case DataType::I32:
    scatterColumn<Tile, std::int32_t>(entry, permutation, first, count, src);
    return;
  case DataType::I64:
    scatterColumn<Tile, std::int64_t>(entry, permutation, first, count, src);
    return;
  }
  logError() << "expr: column" << entry.name << "has an unhandled storage type.";
}

std::int32_t readGroup(const DataEntry& entry, std::size_t index) {
  switch (entry.datatype) {
  case DataType::I32:
    return entry.getValue<std::int32_t>(index);
  case DataType::I64:
    return static_cast<std::int32_t>(entry.getValue<std::int64_t>(index));
  case DataType::F32:
    return static_cast<std::int32_t>(entry.getValue<float>(index));
  case DataType::F64:
    return static_cast<std::int32_t>(entry.getValue<double>(index));
  }
  return 0;
}

} // namespace

Binding Binding::bind(const Program& program, const DataTable& table) {
  Binding binding;
  binding.numPoints_ = table.numPoints();
  binding.computeType_ = program.computeType();

  if (binding.numPoints_ == 0) {
    throw std::invalid_argument("expr: cannot bind a program to an empty point set");
  }

  // Name -> column. A duplicate name is rejected rather than resolved to the
  // first or the last match: which one wins is invisible at the call site, and
  // "the name silently resolved to the wrong column" is the defect class this
  // whole layer is built to remove.
  std::unordered_map<std::string, std::size_t> byName;
  for (std::size_t i = 0; i < table.dataEntries().size(); ++i) {
    const auto& entry = table.dataEntries()[i];
    if (!byName.emplace(entry.name, i).second) {
      throw std::invalid_argument("expr: the data table declares '" + entry.name + "' twice");
    }
  }

  // A channel that is both an input and a state has two different meanings in
  // the lowering, which resolves state first and would silently ignore the
  // column. Caught here because Binding is where the column actually exists.
  for (const auto& state : program.state()) {
    for (const auto& input : program.inputs()) {
      if (input.name == state.name) {
        throw std::invalid_argument("expr: '" + state.name +
                                    "' is declared both as an input and as a state");
      }
    }
  }

  // Inputs, in Program::inputs() order, because Opcode::LoadInput addresses the
  // gathered tile by exactly that index. The vector is dense: every input must
  // resolve, so slot == position and the backend needs no indirection.
  binding.inputs_.reserve(program.inputs().size());
  for (std::size_t i = 0; i < program.inputs().size(); ++i) {
    const auto& spec = program.inputs()[i];
    const auto found = byName.find(spec.name);
    if (found == byName.end()) {
      throw std::invalid_argument("expr: the program reads '" + spec.name +
                                  "', which the data table does not provide");
    }
    const auto& entry = table.dataEntries()[found->second];
    if (entry.direction == Direction::Out) {
      throw std::invalid_argument("expr: the program reads '" + spec.name +
                                  "', which the data table offers for output only");
    }
    ColumnBinding column;
    column.entry = found->second;
    column.slot = static_cast<int>(i);
    column.tableType = entry.datatype;
    column.computed = entry.setter == nullptr;
    binding.inputs_.push_back(column);
  }

  binding.outputs_.reserve(program.outputs().size());
  for (std::size_t i = 0; i < program.outputs().size(); ++i) {
    const auto& spec = program.outputs()[i];
    const auto found = byName.find(spec.name);
    if (found == byName.end()) {
      throw std::invalid_argument("expr: the program writes '" + spec.name +
                                  "', which the data table does not provide");
    }
    const auto& entry = table.dataEntries()[found->second];
    if (entry.direction == Direction::In) {
      throw std::invalid_argument("expr: the program writes '" + spec.name +
                                  "', which the data table offers for input only");
    }
    // Direction and writability are two different questions. A column bound
    // through bindViewConst or bindMemberViewConst can carry Direction::InOut
    // and still have no setter, in which case the scatter would call a null
    // std::function at the first tile. Checked here, once, rather than found at
    // run time on one rank.
    if (entry.setter == nullptr) {
      throw std::invalid_argument("expr: the program writes '" + spec.name +
                                  "', which is bound as a read-only view");
    }
    ColumnBinding column;
    column.entry = found->second;
    column.slot = static_cast<int>(i);
    column.tableType = entry.datatype;
    column.computed = false;
    binding.outputs_.push_back(column);
  }

  binding.buildGroupRanges(program, table);
  return binding;
}

void Binding::buildGroupRanges(const Program& program, const DataTable& table) {
  // Only a program that actually reads `group` gets partitioned. Sorting a
  // point set the kernel will not branch on would cost a permutation indirection
  // in every gather for nothing.
  std::size_t groupEntry = 0;
  bool hasGroup = false;
  for (std::size_t i = 0; i < program.inputs().size(); ++i) {
    if (program.inputs()[i].name == GroupChannelName) {
      groupEntry = inputs_[i].entry;
      hasGroup = true;
      break;
    }
  }
  if (!hasGroup) {
    return;
  }

  const DataEntry& entry = table.dataEntries()[groupEntry];
  std::vector<std::int32_t> groups(numPoints_);
  for (std::size_t p = 0; p < numPoints_; ++p) {
    groups[p] = readGroup(entry, p);
  }

  std::vector<std::size_t> order(numPoints_);
  std::iota(order.begin(), order.end(), std::size_t{0});
  // Stable, so the permutation is a function of the data and not of the sort's
  // internal choices -- two ranks with the same points must agree, or the state
  // slots they carry across calls stop describing the same points.
  std::stable_sort(order.begin(), order.end(), [&groups](std::size_t a, std::size_t b) {
    return groups[a] < groups[b];
  });

  groupRanges_.clear();
  std::size_t begin = 0;
  for (std::size_t k = 1; k <= numPoints_; ++k) {
    if (k == numPoints_ || groups[order[k]] != groups[order[begin]]) {
      groupRanges_.push_back(GroupRange{groups[order[begin]], begin, k});
      begin = k;
    }
  }

  // The common cases -- one group, or a point set already grouped by
  // construction -- leave the order untouched. Storing an identity permutation
  // would buy nothing and cost an indirection per lane in every gather, so it
  // is dropped and `permutation()` stays empty, exactly as for a program with
  // no group input.
  const bool identity = std::is_sorted(groups.begin(), groups.end());
  if (!identity) {
    permutation_ = std::move(order);
  }

  logInfo() << "expr: partitioned" << numPoints_ << "points into" << groupRanges_.size()
            << "group ranges" << (permutation_.empty() ? "(already ordered)" : "(reordered)");
}

void Binding::allocatePersistent(const Program& program, std::int32_t slotCount) {
  // Idempotent on shape. makeKernel calls this, so building an interpreter and
  // an RTC kernel over one Binding would otherwise reset the state slots in
  // between -- and the differential check that compares two backends is exactly
  // the case that does it.
  if (slotCount == persistentSlotCount_ && !persistent_.empty()) {
    return;
  }
  persistentSlotCount_ = slotCount;
  if (slotCount <= 0) {
    persistent_.clear();
    persistent_.shrink_to_fit();
    return;
  }

  const std::size_t elements = static_cast<std::size_t>(slotCount) * numPoints_;
  const std::size_t width = computeType_ == ComputeType::F32 ? sizeof(float) : sizeof(double);
  persistent_.assign(elements * width, std::byte{0});

  // Only the state slots have a defined initial value; the hoisted ones are
  // written by the Precompute stage before the first run(), and reading one
  // before that is a lowering bug rather than something to paper over with a
  // default here.
  if (computeType_ == ComputeType::F32) {
    initialiseState<float>(program, persistentF32(), numPoints_);
  } else {
    initialiseState<double>(program, persistentF64(), numPoints_);
  }
}

double* Binding::persistentF64() {
  if (computeType_ != ComputeType::F64) {
    logError() << "expr: the program computes in f32; asked for the f64 persistent buffer.";
  }
  return persistent_.empty() ? nullptr : reinterpret_cast<double*>(persistent_.data());
}

float* Binding::persistentF32() {
  if (computeType_ != ComputeType::F32) {
    logError() << "expr: the program computes in f64; asked for the f32 persistent buffer.";
  }
  return persistent_.empty() ? nullptr : reinterpret_cast<float*>(persistent_.data());
}

void Binding::gather(const DataTable& table,
                     std::size_t first,
                     std::size_t count,
                     double* dst) const {
  for (const auto& column : inputs_) {
    gatherOne<double>(column,
                      table,
                      permutation_,
                      first,
                      count,
                      dst + static_cast<std::size_t>(column.slot) * count);
  }
}

void Binding::gather(const DataTable& table,
                     std::size_t first,
                     std::size_t count,
                     float* dst) const {
  for (const auto& column : inputs_) {
    gatherOne<float>(column,
                     table,
                     permutation_,
                     first,
                     count,
                     dst + static_cast<std::size_t>(column.slot) * count);
  }
}

void Binding::scatter(const DataTable& table,
                      std::size_t first,
                      std::size_t count,
                      const double* src) const {
  for (const auto& column : outputs_) {
    scatterOne<double>(column,
                       table,
                       permutation_,
                       first,
                       count,
                       src + static_cast<std::size_t>(column.slot) * count);
  }
}

void Binding::scatter(const DataTable& table,
                      std::size_t first,
                      std::size_t count,
                      const float* src) const {
  for (const auto& column : outputs_) {
    scatterOne<float>(column,
                      table,
                      permutation_,
                      first,
                      count,
                      src + static_cast<std::size_t>(column.slot) * count);
  }
}

} // namespace seissol::expr

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "Reader/Scripting/CompiledReader.h"

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Program.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"
#include "utils/logger.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace seissol::reader::scripting {

namespace {

std::vector<std::string> namesOf(const std::vector<expr::VarSpec>& specs) {
  std::vector<std::string> names;
  names.reserve(specs.size());
  for (const auto& spec : specs) {
    names.push_back(spec.name);
  }
  return names;
}

double relativeDeviation(double a, double b) {
  if (std::isnan(a) && std::isnan(b)) {
    return 0.0;
  }
  if (a == b) {
    return 0.0;
  }
  const double scale = std::max({std::fabs(a), std::fabs(b), 1.0});
  return std::fabs(a - b) / scale;
}

} // namespace

CompiledReader::CompiledReader(expr::Program program,
                               std::unique_ptr<DataReader> reference,
                               CompiledReaderOptions options,
                               datafield::GridStore* grids)
    : program_(std::move(program)), reference_(std::move(reference)), options_(std::move(options)),
      inputs_(namesOf(program_.inputs())), outputs_(namesOf(program_.outputs())),
      grids_(grids != nullptr ? grids : &datafield::sharedGridStore()) {}

CompiledReader::~CompiledReader() = default;

void CompiledReader::prepare(const DataTable& table) {
  if (preparedFor_ == &table) {
    return;
  }
  if (preparedFor_ != nullptr) {
    // Not silent: a rebind re-runs the precompute stage and resets every state
    // slot to its initial value, so a consumer alternating between two tables
    // loses the history it thinks it is accumulating.
    logWarning() << "expr: rebinding a compiled reader to a different table; state slots reset.";
  }

  binding_ = expr::Binding::bind(program_, table);
  kernel_ = expr::makeKernel(program_, binding_, *grids_, options_.backend);
  if (kernel_ == nullptr) {
    throw std::runtime_error("expr: no usable backend for the compiled reader");
  }
  kernel_->precompute(table);
  preparedFor_ = &table;

  fallback_ = !differentialCheck(table);
}

bool CompiledReader::differentialCheck(const DataTable& table) {
  if (reference_ == nullptr || options_.differentialSamples == 0) {
    return true;
  }
  if (binding_.outputs().empty()) {
    return true;
  }

  // The outputs have to be readable to be compared. Direction::Out columns bound
  // through a view carry an accessor, but nothing in the ABI guarantees it, and
  // a check that silently passes when it cannot read is worse than no check.
  for (const auto& column : binding_.outputs()) {
    if (table.dataEntries()[column.entry].accessor == nullptr) {
      logWarning() << "expr: skipping the differential check --"
                   << table.dataEntries()[column.entry].name.c_str()
                   << "is bound write-only, so the two paths cannot be compared.";
      return true;
    }
  }

  const std::size_t numPoints = table.numPoints();
  const std::size_t stride = std::max<std::size_t>(1, numPoints / options_.differentialSamples);

  reference_->prepare(table);
  reference_->call(table);

  // Snapshot BEFORE the kernel runs, because both write the same columns.
  std::vector<std::vector<double>> expected(binding_.outputs().size());
  for (std::size_t c = 0; c < binding_.outputs().size(); ++c) {
    const auto& entry = table.dataEntries()[binding_.outputs()[c].entry];
    for (std::size_t p = 0; p < numPoints; p += stride) {
      expected[c].push_back(entry.getValueAs<double>(p));
    }
  }

  kernel_->run(table);

  double worst = 0.0;
  std::string worstName;
  std::size_t worstPoint = 0;
  for (std::size_t c = 0; c < binding_.outputs().size(); ++c) {
    const auto& entry = table.dataEntries()[binding_.outputs()[c].entry];
    std::size_t k = 0;
    for (std::size_t p = 0; p < numPoints; p += stride, ++k) {
      const double deviation = relativeDeviation(entry.getValueAs<double>(p), expected[c][k]);
      if (deviation > worst) {
        worst = deviation;
        worstName = entry.name;
        worstPoint = p;
      }
    }
  }

  if (worst > options_.differentialTolerance) {
    logWarning() << "expr: the compiled program disagrees with the interpreted reader by" << worst
                 << "(relative) on" << worstName.c_str() << "at point" << worstPoint
                 << "-- falling back to the interpreted path.";
    // Leave the table holding the values that are trusted, not the ones that
    // are not: this is prepare(), and a consumer is allowed to read afterwards.
    reference_->call(table);
    return false;
  }

  logInfo() << "expr: differential check passed over" << (numPoints / stride + 1) << "of"
            << numPoints << "points; worst relative deviation" << worst << ".";
  return true;
}

void CompiledReader::call(const DataTable& table) {
  if (preparedFor_ != &table) {
    prepare(table);
    // prepare() already evaluated the table on both paths for the check, but
    // running again is cheap next to the binding and keeps call()'s contract
    // simple: after call(), the outputs come from the path in use.
  }
  if (fallback_) {
    reference_->call(table);
    return;
  }
  kernel_->run(table);
}

} // namespace seissol::reader::scripting

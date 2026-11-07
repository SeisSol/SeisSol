// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_DATAFIELDPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_DATAFIELDPARAMETERS_H_

#include "ParameterReader.h"
#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <optional>

namespace seissol::initializer::parameters {

/// Settings for externally sampled data grids (velocity models, time-dependent
/// boundary data). These are defaults: a grid may override the scheme and the
/// out-of-domain behaviour per grid, since a run can legitimately want cubic for
/// a smooth velocity model and nearest for a categorical one in the same file
/// set.
struct DatafieldParameters {
  /// Total memory, per rank, for the resident time windows of all
  /// time-dependent grids. The sampling synchronisation interval follows from
  /// this: more memory buys a longer window buys fewer reloads.
  std::size_t windowMemory{reader::datafield::DefaultWindowMemoryBudget};
  /// Explicit window size in slices, overriding the memory derivation. Absent
  /// means derive.
  std::optional<std::size_t> residentSlices;
  /// Fallback spacing for a time axis whose file does not state one. Zero means
  /// no fallback, and such a file is an error rather than a guess.
  double defaultTimeSpacing{0.0};
  reader::datafield::Interpolation interpolation{reader::datafield::Interpolation::Linear};
  reader::datafield::BoundaryMode boundary{reader::datafield::BoundaryMode::Clamp};
};

DatafieldParameters readDatafieldParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_DATAFIELDPARAMETERS_H_

// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DatafieldParameters.h"

#include "Initializer/Parameters/ParameterReader.h"
#include "Reader/Datafield/Grid.h"

#include <cstddef>
#include <optional>
#include <string>
#include <utils/logger.h>

namespace seissol::initializer::parameters {

namespace {
using seissol::reader::datafield::BoundaryMode;
using seissol::reader::datafield::Interpolation;
} // namespace

DatafieldParameters readDatafieldParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("datafield");

  // Expressed in mebibytes rather than bytes: the value is a human's judgement
  // about how much of the node a velocity model may occupy, and nobody has an
  // opinion accurate to the byte.
  const auto windowMemoryMiB =
      reader->readWithDefault("windowmemory",
                              static_cast<double>(reader::datafield::DefaultWindowMemoryBudget) /
                                  static_cast<double>(1ULL << 20));
  if (!(windowMemoryMiB > 0.0)) {
    logError() << "datafield: windowmemory must be positive, got" << windowMemoryMiB << "MiB.";
  }
  const auto windowMemory =
      static_cast<std::size_t>(windowMemoryMiB * static_cast<double>(1ULL << 20));

  // Zero means "derive from windowmemory". An explicit count is for the cases
  // where someone knows the file and wants a specific number of slices.
  const auto residentSlicesRaw = reader->readWithDefault("residentslices", 0);
  const auto residentSlices = (residentSlicesRaw > 0)
                                  ? std::optional<std::size_t>(residentSlicesRaw)
                                  : std::optional<std::size_t>{};

  const auto defaultTimeSpacing = reader->readWithDefault("defaulttimespacing", 0.0);
  if (defaultTimeSpacing < 0.0) {
    logError() << "datafield: defaulttimespacing must not be negative, got" << defaultTimeSpacing
               << ".";
  }

  const auto interpolation =
      reader->readWithDefaultStringEnum<Interpolation>("interpolation",
                                                       "linear",
                                                       {
                                                           {"nearest", Interpolation::Nearest},
                                                           {"linear", Interpolation::Linear},
                                                           {"cubic", Interpolation::Cubic},
                                                       });
  if (interpolation == Interpolation::Cubic) {
    // Worth saying out loud: cubic convolution overshoots across a sharp
    // material contrast by roughly 7% of the jump, which can put rho or mu
    // outside the sampled range without anything downstream noticing.
    logInfo() << "datafield: cubic interpolation is not bound-preserving. Across a sharp"
              << "contrast it overshoots, so sampled values may fall outside the range present"
              << "in the file.";
  }

  const auto boundary =
      reader->readWithDefaultStringEnum<BoundaryMode>("boundary",
                                                      "clamp",
                                                      {
                                                          {"clamp", BoundaryMode::Clamp},
                                                          {"fail", BoundaryMode::Fail},
                                                      });

  return DatafieldParameters{
      windowMemory, residentSlices, defaultTimeSpacing, interpolation, boundary};
}
} // namespace seissol::initializer::parameters

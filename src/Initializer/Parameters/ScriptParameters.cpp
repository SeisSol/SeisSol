// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ScriptParameters.h"

#include "Expr/Program.h"
#include "Initializer/Parameters/ParameterReader.h"

#include <cstdint>
#include <optional>

namespace seissol::initializer::parameters {

namespace {
enum class Precision : std::uint8_t { Auto, F32, F64 };
} // namespace

ScriptParameters readScriptParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("scripts");
  ScriptParameters parameters;
  const auto precision = reader->readWithDefaultStringEnum<Precision>("precision",
                                                                      "auto",
                                                                      {
                                                                          {"auto", Precision::Auto},
                                                                          {"f32", Precision::F32},
                                                                          {"f64", Precision::F64},
                                                                      });
  if (precision == Precision::F32) {
    parameters.precision = expr::ComputeType::F32;
  } else if (precision == Precision::F64) {
    parameters.precision = expr::ComputeType::F64;
  }
  parameters.strict = reader->readWithDefault("strict", false);
  return parameters;
}

} // namespace seissol::initializer::parameters

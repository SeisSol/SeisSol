// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_SCRIPTPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_SCRIPTPARAMETERS_H_

#include "Expr/Program.h"
#include "ParameterReader.h"

#include <optional>

namespace seissol::initializer::parameters {

/// The section `&scripts`: how a run evaluates the scripts it evaluates per time step or per output
/// (the boundaries, the slip rates of FL = 36, the derived outputs). In the parameter file rather
/// than in the environment, since they change the values a run computes.
struct ScriptParameters {
  /// Their precision (`precision`): none for that of each configuration (`auto`), else binary32
  /// (`f32`) or binary64 (`f64`).
  std::optional<expr::ComputeType> precision;
  /// Whether a device computes them as written (`strict = 1`): without contracting a product and a
  /// sum into a fused multiply-add, and with correctly rounded division and square root through
  /// OpenCL.
  bool strict{false};
};

ScriptParameters readScriptParameters(ParameterReader* baseReader);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_SCRIPTPARAMETERS_H_

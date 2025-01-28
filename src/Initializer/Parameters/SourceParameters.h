// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_SOURCEPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_SOURCEPARAMETERS_H_

#include <string>

#include "ParameterReader.h"

namespace seissol::initializer::parameters {

enum class PointSourceType : int { None = 0, NrfSource = 42, FsrmSource = 50 };

struct SourceParameters {
  PointSourceType type;
  std::string fileName;
};

SourceParameters readSourceParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_SOURCEPARAMETERS_H_

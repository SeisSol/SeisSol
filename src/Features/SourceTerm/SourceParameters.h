// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_SOURCETERM_SOURCEPARAMETERS_H_
#define SEISSOL_SRC_FEATURES_SOURCETERM_SOURCEPARAMETERS_H_

#include <string>

namespace seissol::initializer::parameters {

enum class PointSourceType : int { None = 0, NrfSource = 42, FsrmSource = 50 };

struct SourceParameters {
  PointSourceType type{PointSourceType::None};
  std::string fileName;
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_FEATURES_SOURCETERM_SOURCEPARAMETERS_H_

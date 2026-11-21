// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_TIMEMIRROR_ITMPARAMETERS_H_
#define SEISSOL_SRC_FEATURES_TIMEMIRROR_ITMPARAMETERS_H_

namespace seissol::initializer::parameters {

enum class ReflectionType { BothWaves = 1, BothWavesVelocity, Pwave, Swave };

struct ITMParameters {
  bool itmEnabled{false};
  double itmStartingTime{0.0};
  double itmDuration{0.0};
  double itmVelocityScalingFactor{1.0};
  ReflectionType itmReflectionType{ReflectionType::BothWaves};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_FEATURES_TIMEMIRROR_ITMPARAMETERS_H_

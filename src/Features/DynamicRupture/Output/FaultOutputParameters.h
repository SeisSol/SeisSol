// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_FAULTOUTPUTPARAMETERS_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_FAULTOUTPUTPARAMETERS_H_

#include "Core/IO/WriterParameters.h"
#include "WavePropagation/Output/ReceiverParameters.h"

#include <array>
#include <cstddef>
#include <optional>
#include <string>

namespace seissol::initializer::parameters {

/// The outputs of a fault: none, its on-fault receivers (pickpoints), its whole surface
/// (elementwise), or both. As a setting (outputpointtype), only its values count; the fault output
/// tells its two outputs apart by AtPickpoint and Elementwise.
enum class OutputType : int {
  None = 0,
  AtPickpoint = 3,
  Elementwise = 4,
  AtPickpointAndElementwise = 5
};

enum class FaultRefinement { Triple = 1, Quad = 2, None = 3 };

struct ElementwiseFaultParameters {
  double printTimeIntervalSec{1.0};
  std::array<bool, 12> outputMask{true, true, true, true};
  FaultRefinement refinementStrategy{FaultRefinement::Quad};
  int refinement{2};
  int vtkorder{-1};
  // none if the parameters choose none (see TimeSeriesMode)
  std::optional<TimeSeriesMode> timeSeries;
};

struct PickpointParameters {
  int printTimeInterval{1};
  double writeInterval{VeryLongTime};
  std::array<bool, 12> outputMask{true, true, true};
  std::optional<std::string> pickpointFileName;
  bool collectiveio{false};
  ReceiverOutputFormat format{ReceiverOutputFormat::Csv};
  //! @brief As ReceiverOutputParameters::samplechunk, for the on-fault receivers.
  std::size_t samplechunk{0};
  /**
   * @brief The time between two samples of the on-fault receivers (`pickdt`); zero samples every
   * `printTimeInterval` time steps of the finest time cluster with fault faces instead.
   */
  double samplingInterval{0};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_FAULTOUTPUTPARAMETERS_H_

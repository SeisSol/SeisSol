// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_DRPARAMETERS_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_DRPARAMETERS_H_

#include "Core/Geometry/FaultOrientation.h"
#include "Features/DynamicRupture/Typedefs.h"

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace seissol::initializer::parameters {

enum class OutputType : int {
  None = 0,
  AtPickpoint = 3,
  Elementwise = 4,
  AtPickpointAndElementwise = 5
};

enum class SlipRateOutputType : int {
  VelocityDifference = 0,
  TractionsAndFailure = 1,
};

struct DRParameters {
  bool isDynamicRuptureEnabled{true};
  bool isThermalPressureOn{false};
  bool isFrictionEnergyRequired{false};
  bool isCheckAbortCriteraEnabled{false};
  bool energiesFromAcrossFaultVelocities{false};
  OutputType outputPointType{OutputType::None};
  RefPointMethod refPointMethod{RefPointMethod::Point};
  SlipRateOutputType slipRateOutputType{SlipRateOutputType::TractionsAndFailure};
  seissol::dr::misc::FrictionLawType frictionLawType{seissol::dr::misc::FrictionLawType::NoFault};
  double healingThreshold{-1.0};
  std::array<double, seissol::dr::MaxNucleations> t0{};
  std::array<double, seissol::dr::MaxNucleations> s0{};
  double tpProxyExponent{1.0 / 3.0};
  double rsF0{0.0};
  double rsB{0.0};
  double rsSr0{0.0};
  double rsInitialSlipRate1{0.0};
  double rsInitialSlipRate2{0.0};
  double muW{0.0};
  double thermalDiffusivity{0.0};
  double heatCapacity{0.0};
  double undrainedTPResponse{0.0};
  double initialTemperature{0.0};
  double initialPressure{0.0};
  double vStar{0.0}; // Prakash-Clifton regularization parameter
  double prakashLength{0.0};
  std::string faultFileName;
  // per fused simulation of the configuration of the run
  std::vector<std::optional<std::string>> faultFileNames;
  Eigen::Vector3d referencePoint{Eigen::Vector3d::Zero()};
  double terminatorSlipRateThreshold{std::numeric_limits<double>::max()};
  double etaDamp{1.0};
  double etaDampEnd{std::numeric_limits<double>::infinity()};
  std::uint32_t nucleationCount{1};
  std::uint32_t rsMaxNumberSlipRateUpdates{60};
  std::uint32_t rsNumberStateVariableUpdates{10};
  double rsSlipRateTolerance{1e-8};
  double rsStateTolerance{1e-8};
  /// points on the logarithmic wavenumber grid the thermal pressurization is integrated on
  std::uint32_t tpGridPoints{seissol::dr::misc::DefaultTpGridPoints};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_DRPARAMETERS_H_

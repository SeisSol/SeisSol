// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_DRPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_DRPARAMETERS_H_

#include "DynamicRupture/Typedefs.h"
#include "ParameterReader.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <cstdint>
#include <numeric>
#include <string>

namespace seissol::initializer::parameters {

enum class OutputType : int {
  None = 0,
  AtPickpoint = 3,
  Elementwise = 4,
  AtPickpointAndElementwise = 5
};

enum class RefPointMethod : int { Point = 0, Normal = 1 };

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
  OutputType outputPointType{3};
  RefPointMethod refPointMethod{0};
  SlipRateOutputType slipRateOutputType{1};
  seissol::dr::misc::FrictionLawType frictionLawType{0};
  double healingThreshold{-1.0};
  std::array<double, seissol::dr::MaxNucleations> t0{};
  std::array<double, seissol::dr::MaxNucleations> s0{};
  double tpProxyExponent{0.0};
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
  std::array<std::optional<std::string>, seissol::multisim::NumSimulations> faultFileNames;
  Eigen::Vector3d referencePoint;
  double terminatorSlipRateThreshold{0.0};
  double etaDamp{1.0};
  double etaDampEnd{std::numeric_limits<double>::infinity()};
  std::uint32_t nucleationCount{0};
  std::uint32_t rsMaxNumberSlipRateUpdates{60};
  std::uint32_t rsNumberStateVariableUpdates{10};
  double rsNewtonTolerance{1e-8};
  double rsStateTolerance{1e-8};
};

DRParameters readDRParameters(ParameterReader* baseReader);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_DRPARAMETERS_H_

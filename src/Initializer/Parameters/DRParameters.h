// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_DRPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_DRPARAMETERS_H_

#include <string>

#include <Eigen/Dense>

#include "Kernels/Precision.h"
#include "ParameterReader.h"

#include <Solver/MultipleSimulations.h>

#include <numeric>

namespace seissol::initializer::parameters {

constexpr std::size_t MaxNucleactions = 16;

/**
 * Stores the different types of friction laws
 * The values resemble the identifiers used in the old fortran implementation
 */
enum class FrictionLawType : unsigned int {
  NoFault = 0,
  LinearSlipWeakeningLegacy = 2,
  LinearSlipWeakening = 16,
  LinearSlipWeakeningBimaterial = 6,
  LinearSlipWeakeningTPApprox = 1058,
  RateAndStateAgingLaw = 3,
  RateAndStateSlipLaw = 4,
  RateAndStateFastVelocityWeakening = 103,
  ImposedSlipRatesYoffe = 33,
  ImposedSlipRatesGaussian = 34,
  ImposedSlipRatesDelta = 35,
  RateAndStateSevereVelocityWeakening = 7,
  RateAndStateAgingNucleation = 101,
};

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

// TODO: maybe make Real<Cfg> again?
struct DRParameters {
  bool isDynamicRuptureEnabled{true};
  bool isThermalPressureOn{false};
  bool isFrictionEnergyRequired{false};
  bool isCheckAbortCriteraEnabled{false};
  bool energiesFromAcrossFaultVelocities{false};
  OutputType outputPointType{3};
  RefPointMethod refPointMethod{0};
  SlipRateOutputType slipRateOutputType{1};
  FrictionLawType frictionLawType{0};
  double healingThreshold{-1.0};
  std::array<double, MaxNucleactions> t0{};
  std::array<double, MaxNucleactions> s0{};
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
  double etaHack{1.0};
  double etaStop{std::numeric_limits<double>::infinity()};
  unsigned nucleationCount{0};
};

DRParameters readDRParameters(ParameterReader* baseReader);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_DRPARAMETERS_H_

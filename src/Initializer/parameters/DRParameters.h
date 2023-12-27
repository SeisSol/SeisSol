#ifndef SEISSOL_DR_PARAMETERS_H
#define SEISSOL_DR_PARAMETERS_H

#include <string>

#include <Eigen/Dense>

#include "ParameterReader.h"

namespace seissol::initializers::parameters {

/**
 * Stores the different types of friction laws
 * The values resemble the identifiers used in the old fortran implementation
 */
enum class FrictionLawType : unsigned int {
  NoFault = 0,
  LinearSlipWeakening = 16,
  LinearSlipWeakeningBimaterial = 6,
  RateAndStateAgingLaw = 3,
  RateAndStateSlipLaw = 4,
  RateAndStateFastVelocityWeakening = 103,
  ImposedSlipRatesYoffe = 33,
  ImposedSlipRatesGaussian = 34,
  RateAndStateVelocityWeakening = 7,
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

struct DRParameters {
  bool isDynamicRuptureEnabled{true};
  bool isRfOutputOn{false};
  bool isDsOutputOn{false};
  bool isThermalPressureOn{false};
  bool isFrictionEnergyRequired{false};
  OutputType outputPointType{3};
  RefPointMethod refPointMethod{0};
  SlipRateOutputType slipRateOutputType{1};
  FrictionLawType frictionLawType{0};
  double t0{0.0};
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
  std::string faultFileName{""};
  Eigen::Vector3d referencePoint;
};

inline DRParameters readDRParameters(ParameterReader& baseReader) {
  auto reader = baseReader.readSubNode("dynamicrupture");

  const double xref = reader.readWithDefault("xref", 0.0);
  const double yref = reader.readWithDefault("yref", 0.0);
  const double zref = reader.readWithDefault("zref", 0.0);
  const Eigen::Vector3d referencePoint = {xref, yref, zref};

  const auto refPointMethod = reader.readWithDefaultEnum<RefPointMethod>(
      "refpointmethod", RefPointMethod::Point, {RefPointMethod::Point, RefPointMethod::Normal});

  const auto outputPointType =
      reader.readWithDefaultEnum<OutputType>("outputpointtype",
                                             OutputType::None,
                                             {OutputType::None,
                                              OutputType::AtPickpoint,
                                              OutputType::Elementwise,
                                              OutputType::AtPickpointAndElementwise});
  const auto frictionLawType = reader.readWithDefaultEnum<FrictionLawType>(
      "fl",
      FrictionLawType::NoFault,
      {FrictionLawType::NoFault,
       FrictionLawType::LinearSlipWeakening,
       FrictionLawType::LinearSlipWeakeningBimaterial,
       FrictionLawType::RateAndStateAgingLaw,
       FrictionLawType::RateAndStateSlipLaw,
       FrictionLawType::RateAndStateFastVelocityWeakening,
       FrictionLawType::ImposedSlipRatesYoffe,
       FrictionLawType::ImposedSlipRatesGaussian,
       FrictionLawType::RateAndStateVelocityWeakening,
       FrictionLawType::RateAndStateAgingNucleation});
  auto slipRateOutputType = reader.readWithDefaultEnum<SlipRateOutputType>(
      "sliprateoutputtype",
      SlipRateOutputType::TractionsAndFailure,
      {SlipRateOutputType::VelocityDifference, SlipRateOutputType::TractionsAndFailure});
  if (((frictionLawType == FrictionLawType::ImposedSlipRatesYoffe) or
       (frictionLawType == FrictionLawType::ImposedSlipRatesGaussian)) and
      (slipRateOutputType == SlipRateOutputType::TractionsAndFailure)) {
    logWarning(seissol::MPI::mpi.rank())
        << "SlipRateOutputType=1 is incompatible with imposed slip rates friction laws, "
           "switching to SlipRateOutputType=0";
    slipRateOutputType = SlipRateOutputType::VelocityDifference;
  }
  const auto backgroundType = reader.readWithDefault("backgroundtype", 0);
  const auto isThermalPressureOn = reader.readWithDefault("thermalpress", false);
  const auto t0 = reader.readWithDefault("t_0", 0.0);

  auto readIfRequired = [&reader](const std::string& name, bool required) {
    double value;
    if (required) {
      std::string failString = "Did not find parameter " + name;
      value = reader.readOrFail<double>(name, failString);
    } else {
      reader.markUnused(name);
    }
    return value;
  };

  const bool isRateAndState =
      (frictionLawType == FrictionLawType::RateAndStateAgingLaw) or
      (frictionLawType == FrictionLawType::RateAndStateSlipLaw) or
      (frictionLawType == FrictionLawType::RateAndStateVelocityWeakening) or
      (frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto rsF0 = readIfRequired("rs_f0", isRateAndState);
  const auto rsB = readIfRequired("rs_b", isRateAndState);
  const auto rsSr0 = readIfRequired("rs_sr0", isRateAndState);
  const auto rsInitialSlipRate1 = readIfRequired("rs_inisliprate1", isRateAndState);
  const auto rsInitialSlipRate2 = readIfRequired("rs_inisliprate2", isRateAndState);

  const auto muW = readIfRequired(
      "rs_muw", frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto thermalDiffusivity = readIfRequired("tp_thermaldiffusivity", isThermalPressureOn);
  const auto heatCapacity = readIfRequired("tp_heatcapacity", isThermalPressureOn);
  const auto undrainedTPResponse = readIfRequired("tp_undrainedtpresponse", isThermalPressureOn);
  const auto initialTemperature = readIfRequired("tp_initemp", isThermalPressureOn);
  const auto initialPressure = readIfRequired("tp_inipressure", isThermalPressureOn);

  const bool isBiMaterial = frictionLawType == FrictionLawType::LinearSlipWeakeningBimaterial;
  const auto vStar = readIfRequired("pc_vstar", isBiMaterial);
  const auto prakashLength = readIfRequired("pc_prakashlength", isBiMaterial);

  const std::string faultFileName = reader.readWithDefault("modelfileName", std::string(""));

  // TODO: isDsOn, isRfOn, isFrictionEnergyOn

  // if there is no fileName given for the fault, assume that we do not use dynamic rupture
  const bool isDynamicRuptureEnabled = faultFileName != "";

  return DRParameters{isDynamicRuptureEnabled,
                      true, // TODO
                      true, // TODO
                      isThermalPressureOn,
                      true, // TODO
                      outputPointType,
                      refPointMethod,
                      slipRateOutputType,
                      frictionLawType,
                      t0,
                      rsF0,
                      rsB,
                      rsSr0,
                      rsInitialSlipRate1,
                      rsInitialSlipRate2,
                      muW,
                      thermalDiffusivity,
                      heatCapacity,
                      undrainedTPResponse,
                      initialTemperature,
                      initialPressure,
                      vStar,
                      prakashLength,
                      faultFileName,
                      referencePoint};
}

} // namespace seissol::initializers::parameters
#endif // SEISSOL_PARAMETERS_H

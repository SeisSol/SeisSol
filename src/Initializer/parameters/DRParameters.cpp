#include "DRParameters.h"
#include <cmath>

namespace seissol::initializers::parameters {

DRParameters readDRParameters(ParameterReader& baseReader) {
  auto& reader = baseReader.readSubNode("dynamicrupture");

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
  const auto t0 = static_cast<real>(reader.readWithDefault("t_0", 0.0));

  auto readIfRequired = [&reader](const std::string& name, bool required) {
    real value;
    if (required) {
      std::string failString = "Did not find parameter " + name;
      value = reader.readOrFail<double>(name, failString);
    } else {
      reader.markUnused({name});
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

  const std::string faultFileName = reader.readWithDefault("modelfilename", std::string(""));

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

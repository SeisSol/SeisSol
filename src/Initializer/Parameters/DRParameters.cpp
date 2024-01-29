#include "DRParameters.h"
#include <cmath>

namespace seissol::initializer::parameters {

DRParameters readDRParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("dynamicrupture");

  const double xref = reader->readWithDefault("xref", 0.0);
  const double yref = reader->readWithDefault("yref", 0.0);
  const double zref = reader->readWithDefault("zref", 0.0);
  const Eigen::Vector3d referencePoint = {xref, yref, zref};

  const auto refPointMethod = reader->readWithDefaultEnum<RefPointMethod>(
      "refpointmethod", RefPointMethod::Point, {RefPointMethod::Point, RefPointMethod::Normal});

  const auto outputPointType =
      reader->readWithDefaultEnum<OutputType>("outputpointtype",
                                              OutputType::None,
                                              {OutputType::None,
                                               OutputType::AtPickpoint,
                                               OutputType::Elementwise,
                                               OutputType::AtPickpointAndElementwise});
  const auto frictionLawType = reader->readWithDefaultEnum<FrictionLawType>(
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
  auto slipRateOutputType = reader->readWithDefaultEnum<SlipRateOutputType>(
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
  const auto backgroundType = reader->readWithDefault("backgroundtype", 0);
  const auto isThermalPressureOn = reader->readWithDefault("thermalpress", false);
  const auto t0 = static_cast<real>(reader->readWithDefault("t_0", 0.0));

  const bool isRateAndState =
      (frictionLawType == FrictionLawType::RateAndStateAgingLaw) or
      (frictionLawType == FrictionLawType::RateAndStateSlipLaw) or
      (frictionLawType == FrictionLawType::RateAndStateVelocityWeakening) or
      (frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto rsF0 = reader->readIfRequired<real>("rs_f0", isRateAndState);
  const auto rsB = reader->readIfRequired<real>("rs_b", isRateAndState);
  const auto rsSr0 = reader->readIfRequired<real>("rs_sr0", isRateAndState);
  const auto rsInitialSlipRate1 = reader->readIfRequired<real>("rs_inisliprate1", isRateAndState);
  const auto rsInitialSlipRate2 = reader->readIfRequired<real>("rs_inisliprate2", isRateAndState);

  const auto muW = reader->readIfRequired<real>(
      "rs_muw", frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto thermalDiffusivity =
      reader->readIfRequired<real>("tp_thermaldiffusivity", isThermalPressureOn);
  const auto heatCapacity = reader->readIfRequired<real>("tp_heatcapacity", isThermalPressureOn);
  const auto undrainedTPResponse =
      reader->readIfRequired<real>("tp_undrainedtpresponse", isThermalPressureOn);
  const auto initialTemperature = reader->readIfRequired<real>("tp_initemp", isThermalPressureOn);
  const auto initialPressure = reader->readIfRequired<real>("tp_inipressure", isThermalPressureOn);

  const bool isBiMaterial = frictionLawType == FrictionLawType::LinearSlipWeakeningBimaterial;
  const auto vStar = reader->readIfRequired<real>("pc_vstar", isBiMaterial);
  const auto prakashLength = reader->readIfRequired<real>("pc_prakashlength", isBiMaterial);

  const std::string faultFileName = reader->readWithDefault("modelfilename", std::string(""));

  auto* outputReader = baseReader->readSubNode("output");
  bool isFrictionEnergyRequired = outputReader->readWithDefault("energyoutput", false);

  // if there is no fileName given for the fault, assume that we do not use dynamic rupture
  const bool isDynamicRuptureEnabled = faultFileName != "";

  reader->warnDeprecated({"rf_output_on"});

  return DRParameters{isDynamicRuptureEnabled,
                      isThermalPressureOn,
                      isFrictionEnergyRequired,
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
} // namespace seissol::initializer::parameters

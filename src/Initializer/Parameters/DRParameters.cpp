#include "DRParameters.h"
#include <Initializer/Parameters/ParameterReader.h>
#include <Kernels/Precision.h>
#include <cmath>
#include <limits>
#include <utils/logger.h>

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
       FrictionLawType::LinearSlipWeakeningTPApprox,
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
  const auto isThermalPressureOn = reader->readWithDefault("thermalpress", false);
  const auto healingThreshold =
      static_cast<real>(reader->readWithDefault("lsw_healingthreshold", -1.0));
  const auto t0 = static_cast<real>(reader->readWithDefault("t_0", 0.0));
  const auto tpProxyExponent =
      static_cast<real>(reader->readWithDefault("tpproxyexponent", 1. / 3.));

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

  const auto faultFileName = reader->readPath("modelfilename");

  auto* outputReader = baseReader->readSubNode("output");
  const bool isFrictionEnergyRequired = outputReader->readWithDefault("energyoutput", false);

  auto* abortCriteriaReader = baseReader->readSubNode("abortcriteria");
  const auto terminatorSlipRateThreshold = static_cast<real>(abortCriteriaReader->readWithDefault(
      "terminatorslipratethreshold", std::numeric_limits<real>::infinity()));
  const auto terminatorMaxTimePostRupture = abortCriteriaReader->readWithDefault(
      "terminatormaxtimepostrupture", std::numeric_limits<double>::infinity());
  const bool isCheckAbortCriteraEnabled = std::isfinite(terminatorMaxTimePostRupture);

  // if there is no fileName given for the fault, assume that we do not use dynamic rupture
  const bool isDynamicRuptureEnabled = faultFileName.value_or("") != "";

  const double etaHack = [&]() {
    const auto hackRead1 = reader->read<double>("etahack");
    if (hackRead1.has_value()) {
      return hackRead1.value();
    } else {
      const auto hackRead2 = outputReader->read<double>("etahack");
      if (hackRead2.has_value()) {
        logWarning(seissol::MPI::mpi.rank())
            << "Reading the etahack parameter from the output section is deprecated and may be "
               "removed in a future version of SeisSol. Put the parameter into the dynamicrupture "
               "section instead.";
      }
      return hackRead2.value_or(1.0);
    }
  }();

  reader->warnDeprecated({"rf_output_on", "backgroundtype"});

  return DRParameters{isDynamicRuptureEnabled,
                      isThermalPressureOn,
                      isFrictionEnergyRequired,
                      isCheckAbortCriteraEnabled,
                      outputPointType,
                      refPointMethod,
                      slipRateOutputType,
                      frictionLawType,
                      healingThreshold,
                      t0,
                      tpProxyExponent,
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
                      faultFileName.value_or(""),
                      referencePoint,
                      terminatorSlipRateThreshold,
                      etaHack};
}
} // namespace seissol::initializer::parameters

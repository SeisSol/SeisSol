// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DRParameters.h"
#include <Initializer/Parameters/ParameterReader.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utils/logger.h>
#include <vector>

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
  auto frictionLawType = reader->readWithDefaultEnum<FrictionLawType>(
      "fl",
      FrictionLawType::NoFault,
      {FrictionLawType::NoFault,
       FrictionLawType::LinearSlipWeakening,
       FrictionLawType::LinearSlipWeakeningLegacy,
       FrictionLawType::LinearSlipWeakeningBimaterial,
       FrictionLawType::LinearSlipWeakeningTPApprox,
       FrictionLawType::RateAndStateAgingLaw,
       FrictionLawType::RateAndStateSlipLaw,
       FrictionLawType::RateAndStateFastVelocityWeakening,
       FrictionLawType::ImposedSlipRatesYoffe,
       FrictionLawType::ImposedSlipRatesGaussian,
       FrictionLawType::ImposedSlipRatesDelta,
       FrictionLawType::RateAndStateSevereVelocityWeakening,
       FrictionLawType::RateAndStateAgingNucleation});
  if (frictionLawType == FrictionLawType::LinearSlipWeakeningLegacy) {
    logWarning() << "Using FL=2 for the linear slip weakening friction law is deprecated; consider "
                    "switching it to FL=16";
  }
  auto slipRateOutputType = reader->readWithDefaultEnum<SlipRateOutputType>(
      "sliprateoutputtype",
      SlipRateOutputType::TractionsAndFailure,
      {SlipRateOutputType::VelocityDifference, SlipRateOutputType::TractionsAndFailure});
  if (((frictionLawType == FrictionLawType::ImposedSlipRatesYoffe) or
       (frictionLawType == FrictionLawType::ImposedSlipRatesGaussian) or
       (frictionLawType == FrictionLawType::ImposedSlipRatesDelta)) and
      (slipRateOutputType == SlipRateOutputType::TractionsAndFailure)) {
    logWarning() << "SlipRateOutputType=1 is incompatible with imposed slip rates friction laws, "
                    "switching to SlipRateOutputType=0";
    slipRateOutputType = SlipRateOutputType::VelocityDifference;
  }
  const auto isThermalPressureOn = reader->readWithDefault("thermalpress", false);
  const auto healingThreshold =
      static_cast<double>(reader->readWithDefault("lsw_healingthreshold", -1.0));
  const auto nucleationCount = reader->readWithDefault("nucleationcount", 1U);
  if (nucleationCount > MaxNucleactions) {
    logError() << "You requested more nucleations than supported by this build of SeisSol. Either "
                  "adjust that yourself, or complain to the developers. :)";
  }
  std::array<double, MaxNucleactions> t0;
  std::array<double, MaxNucleactions> s0;
  for (std::size_t i = 0; i < nucleationCount; ++i) {
    const std::string t0name = i == 0 ? "t_0" : ("t" + std::to_string(i + 1) + "_0");
    t0[i] = static_cast<double>(reader->readWithDefault(t0name, 0.0));
    const std::string s0name = i == 0 ? "s_0" : ("s" + std::to_string(i + 1) + "_0");
    s0[i] = static_cast<double>(reader->readWithDefault(s0name, 0.0));
  }
  const auto tpProxyExponent =
      static_cast<double>(reader->readWithDefault("tpproxyexponent", 1. / 3.));

  const bool isRateAndState =
      (frictionLawType == FrictionLawType::RateAndStateAgingLaw) or
      (frictionLawType == FrictionLawType::RateAndStateSlipLaw) or
      (frictionLawType == FrictionLawType::RateAndStateSevereVelocityWeakening) or
      (frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto rsF0 = reader->readIfRequired<double>("rs_f0", isRateAndState);
  const auto rsB = reader->readIfRequired<double>("rs_b", isRateAndState);
  const auto rsSr0 = reader->readIfRequired<double>("rs_sr0", isRateAndState);
  const auto rsInitialSlipRate1 = reader->readIfRequired<double>("rs_inisliprate1", isRateAndState);
  const auto rsInitialSlipRate2 = reader->readIfRequired<double>("rs_inisliprate2", isRateAndState);

  const auto muW = reader->readIfRequiredAlternatives<double>(
      {"rs_muw", "mu_w"}, frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto thermalDiffusivity = reader->readIfRequiredAlternatives<double>(
      {"tp_thermaldiffusivity", "alpha_th"}, isThermalPressureOn);
  const auto heatCapacity =
      reader->readIfRequiredAlternatives<double>({"tp_heatcapacity", "rho_c"}, isThermalPressureOn);
  const auto undrainedTPResponse = reader->readIfRequiredAlternatives<double>(
      {"tp_undrainedtpresponse", "tp_lambda"}, isThermalPressureOn);
  const auto initialTemperature =
      reader->readIfRequiredAlternatives<double>({"tp_initemp", "initemp"}, isThermalPressureOn);
  const auto initialPressure = reader->readIfRequiredAlternatives<double>(
      {"tp_inipressure", "inipressure"}, isThermalPressureOn);

  const bool isBiMaterial = frictionLawType == FrictionLawType::LinearSlipWeakeningBimaterial;
  const auto vStar =
      reader->readIfRequiredAlternatives<double>({"pc_vstar", "v_star"}, isBiMaterial);
  const auto prakashLength =
      reader->readIfRequiredAlternatives<double>({"pc_prakashlength", "l"}, isBiMaterial);

  const auto faultFileName = reader->readPath("modelfilename");

  std::vector<std::optional<std::string>> faultFileNames;

  bool isDynamicRuptureEnabled = false;

  if (!faultFileName.value_or("").empty()) {
    faultFileNames.emplace_back(faultFileName.value());
    isDynamicRuptureEnabled = true;
  }

  for (std::size_t i = 0; i < 64; ++i) {
    const auto fieldname = "modelfilename" + std::to_string(i);
    if (reader->hasField(fieldname)) {
      for (std::size_t j = faultFileNames.size(); j < i; ++j) {
        faultFileNames.emplace_back();
      }
      faultFileNames[i] = reader->read<std::string>(fieldname);
      isDynamicRuptureEnabled = true;
    }
  }

  auto* outputReader = baseReader->readSubNode("output");
  const bool isFrictionEnergyRequired = outputReader->readWithDefault("energyoutput", false);
  const bool energiesFromAcrossFaultVelocities =
      outputReader->readWithDefault("faultenergiesfromacrossfaultvelocities", false);

  auto* abortCriteriaReader = baseReader->readSubNode("abortcriteria");
  const auto terminatorSlipRateThreshold = static_cast<double>(abortCriteriaReader->readWithDefault(
      "terminatorslipratethreshold", std::numeric_limits<double>::infinity()));
  const auto terminatorMaxTimePostRupture = abortCriteriaReader->readWithDefault(
      "terminatormaxtimepostrupture", std::numeric_limits<double>::infinity());
  const bool isCheckAbortCriteraEnabled = std::isfinite(terminatorMaxTimePostRupture);

  const double etaHack = [&]() {
    const auto hackRead1 = reader->read<double>("etahack");
    if (hackRead1.has_value()) {
      return hackRead1.value();
    } else {
      const auto hackRead2 = outputReader->read<double>("etahack");
      if (hackRead2.has_value()) {
        logWarning()
            << "Reading the etahack parameter from the output section is deprecated and may be "
               "removed in a future version of SeisSol. Put the parameter into the dynamicrupture "
               "section instead.";
      }
      return hackRead2.value_or(1.0);
    }
  }();

  const auto hackStop =
      reader->read<double>("etastop").value_or(std::numeric_limits<double>::infinity());

  reader->warnDeprecated({"rf_output_on", "backgroundtype"});

  return DRParameters{isDynamicRuptureEnabled,
                      isThermalPressureOn,
                      isFrictionEnergyRequired,
                      isCheckAbortCriteraEnabled,
                      energiesFromAcrossFaultVelocities,
                      outputPointType,
                      refPointMethod,
                      slipRateOutputType,
                      frictionLawType,
                      healingThreshold,
                      t0,
                      s0,
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
                      faultFileNames,
                      referencePoint,
                      terminatorSlipRateThreshold,
                      etaHack,
                      hackStop,
                      nucleationCount};
}
} // namespace seissol::initializer::parameters

// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DRParameters.h"
#include <Initializer/Parameters/ParameterReader.h>
#include <Kernels/Precision.h>
#include <Solver/MultipleSimulations.h>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
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
      static_cast<real>(reader->readWithDefault("lsw_healingthreshold", -1.0));
  const auto nucleationCount = reader->readWithDefault("nucleationcount", 1U);
  if (nucleationCount > MaxNucleactions) {
    logError() << "You requested more nucleations than supported by this build of SeisSol. Either "
                  "adjust that yourself, or complain to the developers. :)";
  }
  std::array<real, MaxNucleactions> t0;
  std::array<real, MaxNucleactions> s0;
  for (std::size_t i = 0; i < nucleationCount; ++i) {
    const std::string t0name = i == 0 ? "t_0" : ("t" + std::to_string(i + 1) + "_0");
    t0[i] = static_cast<real>(reader->readWithDefault(t0name, 0.0));
    const std::string s0name = i == 0 ? "s_0" : ("s" + std::to_string(i + 1) + "_0");
    s0[i] = static_cast<real>(reader->readWithDefault(s0name, 0.0));
  }
  const auto tpProxyExponent =
      static_cast<real>(reader->readWithDefault("tpproxyexponent", 1. / 3.));

  const bool isRateAndState =
      (frictionLawType == FrictionLawType::RateAndStateAgingLaw) or
      (frictionLawType == FrictionLawType::RateAndStateSlipLaw) or
      (frictionLawType == FrictionLawType::RateAndStateSevereVelocityWeakening) or
      (frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto rsF0 = reader->readIfRequired<real>("rs_f0", isRateAndState);
  const auto rsB = reader->readIfRequired<real>("rs_b", isRateAndState);
  const auto rsSr0 = reader->readIfRequired<real>("rs_sr0", isRateAndState);
  const auto rsInitialSlipRate1 = reader->readIfRequired<real>("rs_inisliprate1", isRateAndState);
  const auto rsInitialSlipRate2 = reader->readIfRequired<real>("rs_inisliprate2", isRateAndState);

  const auto muW = reader->readIfRequiredAlternatives<real>(
      {"rs_muw", "mu_w"}, frictionLawType == FrictionLawType::RateAndStateFastVelocityWeakening);

  const auto thermalDiffusivity = reader->readIfRequiredAlternatives<real>(
      {"tp_thermaldiffusivity", "alpha_th"}, isThermalPressureOn);
  const auto heatCapacity =
      reader->readIfRequiredAlternatives<real>({"tp_heatcapacity", "rho_c"}, isThermalPressureOn);
  const auto undrainedTPResponse = reader->readIfRequiredAlternatives<real>(
      {"tp_undrainedtpresponse", "tp_lambda"}, isThermalPressureOn);
  const auto initialTemperature =
      reader->readIfRequiredAlternatives<real>({"tp_initemp", "initemp"}, isThermalPressureOn);
  const auto initialPressure = reader->readIfRequiredAlternatives<real>(
      {"tp_inipressure", "inipressure"}, isThermalPressureOn);

  const bool isBiMaterial = frictionLawType == FrictionLawType::LinearSlipWeakeningBimaterial;
  const auto vStar = reader->readIfRequiredAlternatives<real>({"pc_vstar", "v_star"}, isBiMaterial);
  const auto prakashLength =
      reader->readIfRequiredAlternatives<real>({"pc_prakashlength", "l"}, isBiMaterial);

  const auto faultFileName = reader->readPath("modelfilename");

  std::array<std::optional<std::string>, seissol::multisim::NumSimulations> faultFileNames;

  bool isDynamicRuptureEnabled = false;

  if (!faultFileName.value_or("").empty()) {
    faultFileNames[0] = faultFileName.value();
    isDynamicRuptureEnabled = true;
  }

  for (std::size_t i = 0; i < faultFileNames.size(); ++i) {
    const auto fieldname = "modelfilename" + std::to_string(i);
    if (reader->hasField(fieldname)) {
      faultFileNames[i] = reader->read<std::string>(fieldname);
      isDynamicRuptureEnabled = true;
    }
  }

  auto* outputReader = baseReader->readSubNode("output");
  const bool isFrictionEnergyRequired = outputReader->readWithDefault("energyoutput", false);
  const bool energiesFromAcrossFaultVelocities =
      outputReader->readWithDefault("faultenergiesfromacrossfaultvelocities", false);

  auto* abortCriteriaReader = baseReader->readSubNode("abortcriteria");
  const auto terminatorSlipRateThreshold = static_cast<real>(abortCriteriaReader->readWithDefault(
      "terminatorslipratethreshold", std::numeric_limits<real>::infinity()));
  const auto terminatorMaxTimePostRupture = abortCriteriaReader->readWithDefault(
      "terminatormaxtimepostrupture", std::numeric_limits<double>::infinity());
  const bool isCheckAbortCriteraEnabled = std::isfinite(terminatorMaxTimePostRupture);

  const double etaDamp = [&]() {
    const auto dampRead1 = reader->read<double>("etadamp");
    if (dampRead1.has_value()) {
      return dampRead1.value();
    } else {
      const auto dampRead2 = reader->read<double>("etahack");
      const auto dampRead3 = outputReader->read<double>("etahack");
      if (dampRead2.has_value() || dampRead3.has_value()) {
        logWarning()
            << "The name \"etahack\" is deprecated and may be "
               "removed in a future version of SeisSol. Put the parameter into the dynamicrupture "
               "section instead and name it \"etadamp\" instead.";
      }
      return dampRead2.value_or(dampRead3.value_or(1.0));
    }
  }();

  const auto etaDampEnd = reader->read<double>("etadampend")
                              .value_or(reader->read<double>("etastop").value_or(
                                  std::numeric_limits<double>::infinity()));

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
                      etaDamp,
                      etaDampEnd,
                      nucleationCount};
}
} // namespace seissol::initializer::parameters

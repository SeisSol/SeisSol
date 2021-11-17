#ifndef SEISSOL_DROUTOUT_DRINITIALIZER_HPP
#define SEISSOL_DROUTOUT_DRINITIALIZER_HPP

#include "DataTypes.hpp"
#include <yaml-cpp/yaml.h>

namespace seissol::dr::output {
class ParametersInitializer {
  public:
  explicit ParametersInitializer(const YAML::Node& userData) : data(userData) {}

  GeneralParamsT getDrGeneralParams() {
    using namespace initializers;
    GeneralParamsT params{};

    if (!data["dynamicrupture"]) {
      throw std::runtime_error("dynamic rupture params. is not provided");
    }

    using namespace seissol::initializers;
    const YAML::Node& drSettings = data["dynamicrupture"];
    auto outputPointID = static_cast<int>(OutputType::AtPickpoint);
    updateIfExists(drSettings, "outputpointtype", outputPointID);
    params.outputPointType = static_cast<OutputType>(outputPointID);

    updateIfExists(drSettings, "sliprateoutputtype", params.slipRateOutputType);
    updateIfExists(drSettings, "fl", params.frictionLawType);
    updateIfExists(drSettings, "backgroundtype", params.backgroundType);
    updateIfExists(drSettings, "rf_output_on", params.isRfOutputOn);
    updateIfExists(drSettings, "ds_output_on", params.isDsOutputOn);
    updateIfExists(drSettings, "magnitude_output_on", params.isMagnitudeOutputOn);
    updateIfExists(drSettings, "energy_rate_output_on", params.isEnergyRateOutputOn);
    updateIfExists(drSettings, "gpwise", params.isGpWiseOutput);
    updateIfExists(drSettings, "thermalpress", params.isTermalPressureOn);
    updateIfExists(drSettings, "energy_rate_printtimeinterval", params.backgroundType);

    using namespace seissol::initializers;
    const YAML::Node& outputParams = data["output"];
    updateIfExists(outputParams, "faultoutputflag", params.faultOutputFlag);
    updateIfExists(outputParams, "outputfile", params.outputFilePrefix);
    updateIfExists(outputParams, "checkpointbackend", params.checkPointBackend);

    const YAML::Node& abortCriteriaParams = data["abortcriteria"];
    updateIfExists(abortCriteriaParams, "endtime", params.endTime);
    updateIfExists(abortCriteriaParams, "maxiteration", params.maxIteration);

#ifdef USE_HDF
    updateIfExists(outputParams, "xdmfwriterbackend", params.xdmfWriterBackend);
#else
    params.xdmfWriterBackend = std::string("posix");
    updateIfExists(outputParams, "xdmfwriterbackend", params.xdmfWriterBackend);
#endif

    return params;
  }

  PickpointParamsT getPickPointParams() {
    using namespace initializers;
    PickpointParamsT ppParams{};

    if (!data["pickpoint"]) {
      throw std::runtime_error("pickpoint output parameters for dynamic rupture is not provided");
    }

    using namespace seissol::initializers;
    const YAML::Node& ppData = data["pickpoint"];
    updateIfExists(ppData, "printtimeinterval", ppParams.printTimeInterval);
    updateIfExists(ppData, "noutpoints", ppParams.numOutputPoints);
    updateIfExists(ppData, "ppfilename", ppParams.ppFileName);

    if (ppData["outputmask"]) {
      convertStringToMask(ppData["outputmask"].as<std::string>(), ppParams.outputMask);
    }

    return ppParams;
  }

  ElementwiseFaultParamsT getElementwiseFaultParams() {
    using namespace initializers;

    ElementwiseFaultParamsT ewParams{};

    if (!data["elementwise"]) {
      throw std::runtime_error(
          "elementwise fault output parameters for dynamic rupture is not provided");
    }

    using namespace seissol::initializers;
    const YAML::Node& ewData = data["elementwise"];
    updateIfExists(ewData, "printtimeinterval", ewParams.printTimeInterval);
    updateIfExists(ewData, "printtimeinterval_sec", ewParams.printTimeIntervalSec);
    updateIfExists(ewData, "printintervalcriterion", ewParams.printIntervalCriterion);
    updateIfExists(ewData, "maxpickstore", ewParams.maxPickStore);
    updateIfExists(ewData, "refinement_strategy", ewParams.refinementStrategy);
    updateIfExists(ewData, "refinement", ewParams.refinement);

    if (ewData["outputmask"]) {
      convertStringToMask(ewData["outputmask"].as<std::string>(), ewParams.outputMask);
    }

    return ewParams;
  }

  private:
  const YAML::Node& data;
};
} // namespace seissol::dr::output
#endif // SEISSOL_DRINITIALIZER_HPP

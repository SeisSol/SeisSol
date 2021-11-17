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
    params.outputPointType =
        static_cast<OutputType>(getParamIfExists(drSettings, "outputpointtype", 3));
    params.slipRateOutputType = getParamIfExists(drSettings, "sliprateoutputtype", 1);
    params.frictionLawType = getParamIfExists(drSettings, "fl", 0);
    params.backgroundType = getParamIfExists(drSettings, "backgroundtype", 0);
    params.isRfOutputOn = getParamIfExists(drSettings, "rf_output_on", false);
    params.isDsOutputOn = getParamIfExists(drSettings, "ds_output_on", false);
    params.isMagnitudeOutputOn = getParamIfExists(drSettings, "magnitude_output_on", false);
    params.isEnergyRateOutputOn = getParamIfExists(drSettings, "energy_rate_output_on", false);
    params.isGpWiseOutput = getParamIfExists(drSettings, "gpwise", false);
    params.isTermalPressureOn = getParamIfExists(drSettings, "thermalpress", false);
    params.backgroundType = getParamIfExists(drSettings, "energy_rate_printtimeinterval", 1);

    using namespace seissol::initializers;
    const YAML::Node& outputParams = data["output"];
    params.faultOutputFlag = getParamIfExists(outputParams, "faultoutputflag", false);
    params.outputFilePrefix = getParamIfExists(outputParams, "outputfile", std::string("data"));
    params.checkPointBackend =
        getParamIfExists(outputParams, "checkpointbackend", std::string("none"));

    const YAML::Node& abortCriteriaParams = data["abortcriteria"];
    params.endTime = getParamIfExists(abortCriteriaParams, "endtime", 15.0);
    params.maxIteration = getParamIfExists(abortCriteriaParams, "maxiteration", 10000000);

#ifdef USE_HDF
    params.xdmfWriterBackend =
        getParamIfExists(outputParams, "xdmfwriterbackend", std::string("hdf5"));
#else
    params.xdmfWriterBackend =
        getParamIfExists(outputParams, "xdmfwriterbackend", std::string("posix"));
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
    ppParams.printTimeInterval = getParamIfExists(ppData, "printtimeinterval", 1);
    ppParams.numOutputPoints = getParamIfExists(ppData, "noutpoints", 0);
    ppParams.ppFileName = getParamIfExists(ppData, "ppfilename", std::string());

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
    ewParams.printTimeInterval = getParamIfExists(ewData, "printtimeinterval", 2);
    ewParams.printTimeIntervalSec = getParamIfExists(ewData, "printtimeinterval_sec", 1.0);
    ewParams.printIntervalCriterion = getParamIfExists(ewData, "printintervalcriterion", 1);
    ewParams.maxPickStore = getParamIfExists(ewData, "maxpickstore", 50);
    ewParams.refinementStrategy = getParamIfExists(ewData, "refinement_strategy", 2);
    ewParams.refinement = getParamIfExists(ewData, "refinement", 2);

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
